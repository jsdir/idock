#include <iostream>
#include <cmath>
#include <cassert>
#include <boost/filesystem/fstream.hpp>
#include "scoring_function.hpp"
#include "array.hpp"
#include "receptor.hpp"

const double receptor::Default_Partition_Granularity = static_cast<double>(3);
const double receptor::Default_Partition_Granularity_Inverse = 1 / Default_Partition_Granularity;

receptor::receptor(const path& p, const array<double, 3>& center, const array<double, 3>& span_, const double grid_granularity) : center(center), grid_granularity(grid_granularity), grid_granularity_inverse(1 / grid_granularity), grid_size({ grid_granularity, grid_granularity, grid_granularity }), grid_size_inverse({ grid_granularity_inverse, grid_granularity_inverse, grid_granularity_inverse }), grid_maps(scoring_function::n), num_probes_product(1)
{
	// The loop may be unrolled by enabling compiler optimization.
	for (size_t i = 0; i < 3; ++i)
	{
		// Slightly expand the user-input span to the nearest multiple of granularity.
		num_grids[i] = static_cast<size_t>(ceil(span_[i] * grid_size_inverse[i]));
		assert(num_grids[i] > 0);
		span[i] = grid_size[i] * num_grids[i];
		num_probes[i] = num_grids[i] + 1;
		num_probes_product *= num_probes[i];

		// Determine the two extreme corners.
		corner0[i] = center[i] - span[i] * static_cast<double>(0.5);
		corner1[i] = corner0[i] + span[i];
		assert(corner0[i] < corner1[i]);

		// Determine the number of partitions.
		num_partitions[i] = static_cast<size_t>(span[i] * Default_Partition_Granularity_Inverse);
		assert(num_partitions[i] > 0);
		partition_size[i] = span[i] / num_partitions[i];
		partition_size_inverse[i] = 1 / partition_size[i];
	}
	partitions.resize(num_partitions[0] * num_partitions[1] * num_partitions[2]);

	// Initialize necessary variables for constructing a receptor.
	atoms.reserve(5000); // A receptor typically consists of <= 5,000 atoms.

	// Initialize helper variables for parsing.
	string residue = "XXXX"; // Current residue sequence, used to track residue change, initialized to a dummy value.
	size_t residue_start; // The starting atom of the current residue.
	string line;

	 // Start parsing.
	for (boost::filesystem::ifstream ifs(p); getline(ifs, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "ATOM  " || record == "HETATM")
		{
			// Parse the residue sequence located at 1-based [23, 26].
			if ((line[25] != residue[3]) || (line[24] != residue[2]) || (line[23] != residue[1]) || (line[22] != residue[0])) // This line is the start of a new residue.
			{
				residue[3] = line[25];
				residue[2] = line[24];
				residue[1] = line[23];
				residue[0] = line[22];
				residue_start = atoms.size();
			}

			// Parse the line.
			atom a(line);

			// Skip unsupported atom types.
			if (a.ad_unsupported()) continue;

			// Skip non-polar hydrogens.
			if (a.is_nonpolar_hydrogen()) continue;

			// For a polar hydrogen, the bonded hetero atom must be a hydrogen bond donor.
			if (a.is_polar_hydrogen())
			{
				for (size_t i = atoms.size(); i > residue_start;)
				{
					atom& b = atoms[--i];
					if (b.is_hetero() && b.is_neighbor(a))
					{
						b.donorize();
						break;
					}
				}
				continue;
			}

			// For a hetero atom, its connected carbon atoms are no longer hydrophobic.
			if (a.is_hetero())
			{
				for (size_t i = atoms.size(); i > residue_start;)
				{
					atom& b = atoms[--i];
					if (!b.is_hetero() && b.is_neighbor(a))
					{
						b.dehydrophobicize();
					}
				}
			}
			// For a carbon atom, it is no longer hydrophobic when connected to a hetero atom.
			else
			{
				for (size_t i = atoms.size(); i > residue_start;)
				{
					const atom& b = atoms[--i];
					if (b.is_hetero() && b.is_neighbor(a))
					{
						a.dehydrophobicize();
						break;
					}
				}
			}

			// Save the atom if and only if its distance to its projection point on the box is within cutoff.
			double r2 = 0;
			for (size_t i = 0; i < 3; ++i)
			{
				if (a.coordinate[i] < corner0[i])
				{
					const double d = a.coordinate[i] - corner0[i];
					r2 += d * d;
				}
				else if (a.coordinate[i] > corner1[i])
				{
					const double d = a.coordinate[i] - corner1[i];
					r2 += d * d;
				}
			}
			if (r2 < scoring_function::cutoff_sqr)
			{
				atoms.push_back(move(a));
			}
		}
		else if (record == "TER   ")
		{
			residue = "XXXX";
		}
	}

	// Allocate each nearby receptor atom to its corresponding partition.
	for (size_t x = 0; x < num_partitions[0]; ++x)
	for (size_t y = 0; y < num_partitions[1]; ++y)
	for (size_t z = 0; z < num_partitions[2]; ++z)
	{
		vector<size_t>& par = partitions[partition_index(array<size_t, 3>{x, y, z})];
		par.reserve(atoms.size());
		const array<size_t, 3> index0 = {{ x,     y,     z     }};
		const array<size_t, 3> index1 = {{ x + 1, y + 1, z + 1 }};
		const array<double, 3> corner0 = partition_corner0(index0);
		const array<double, 3> corner1 = partition_corner0(index1);
		for (size_t i = 0; i < atoms.size(); ++i)
		{
			const atom& a = atoms[i];
			const double proj_dist_sqr = project_distance_sqr(corner0, corner1, a.coordinate);
			if (proj_dist_sqr < scoring_function::cutoff_sqr)
			{
				par.push_back(i);
			}
		}
	}
}

bool receptor::within(const array<double, 3>& coordinate) const
{
	for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
	{
		// Half-open-half-close box, i.e. [corner0, corner1)
		if (coordinate[i] < corner0[i] || corner1[i] <= coordinate[i])
			return false;
	}
	return true;
}

double receptor::project_distance_sqr(const array<double, 3>& corner0, const array<double, 3>& corner1, const array<double, 3>& coordinate) const
{
	// Calculate the projection point of the given coordinate onto the surface of the given box.
	array<double, 3> projection = coordinate; // The loop may be unrolled by enabling compiler optimization.
	for (size_t i = 0; i < 3; ++i)
	{
		if (projection[i] < corner0[i]) projection[i] = corner0[i];
		if (projection[i] > corner1[i]) projection[i] = corner1[i];
	}

	// Check if the distance between the projection and the given coordinate is within cutoff.
	return distance_sqr(projection, coordinate);
}

double receptor::project_distance_sqr(const array<double, 3>& coordinate) const
{
	return project_distance_sqr(corner0, corner1, coordinate);
}

array<double, 3> receptor::grid_corner0(const array<size_t, 3>& index) const
{
	return corner0 + (grid_size * index);
}

array<double, 3> receptor::partition_corner0(const array<size_t, 3>& index) const
{
	return corner0 + (partition_size * index);
}

array<size_t, 3> receptor::grid_index(const array<double, 3>& coordinate) const
{
	array<size_t, 3> index;
	for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
	{
		index[i] = static_cast<size_t>((coordinate[i] - corner0[i]) * grid_size_inverse[i]);
		// Boundary checking is not necessary because the given coordinate is a ligand atom,
		// which has been restricted within the half-open-half-close box [corner0, corner1).
		//if (index[i] == num_grids[i]) index[i] = num_grids[i] - 1;
	}
	return index;
}

array<size_t, 3> receptor::partition_index(const array<double, 3>& coordinate) const
{
	array<size_t, 3> index;
	for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
	{
		index[i] = static_cast<size_t>((coordinate[i] - corner0[i]) * partition_size_inverse[i]);
		// The following condition occurs if and only if coordinate[i] is exactly at the right boundary of the box.
		// In such case, merge it into the last partition.
		// Boundary checking is necessary because the given coordinate is a probe atom.
		if (index[i] == num_partitions[i]) index[i] = num_partitions[i] - 1;
	}
	return index;
}

size_t receptor::grid_index(const array<size_t, 3>& a) const
{
	return num_probes[2] * (num_probes[1] * a[0] + a[1]) + a[2];
}

size_t receptor::partition_index(const array<size_t, 3>& a) const
{
	return num_partitions[2] * (num_partitions[1] * a[0] + a[1]) + a[2];
}

void receptor::populate(const vector<size_t>& xs, const size_t x, const scoring_function& sf)
{
	const size_t num_xs = xs.size();
	vector<double> e(num_xs);

	// For each probe atom of the given X dimension value.
	const size_t num_y_probes = num_probes[1];
	const size_t num_z_probes = num_probes[2];
	for (size_t y = 0; y < num_y_probes; ++y)
	for (size_t z = 0; z < num_z_probes; ++z)
	{
		// Find the possibly interacting receptor atoms via partitions.
		const array<size_t, 3> grid_idx = { { x, y, z } };
		const array<double, 3> probe_coords = grid_corner0(grid_idx);
		const vector<size_t>& receptor_atoms = partitions[partition_index(partition_index(probe_coords))];

		// Accumulate individual free energies for each atom types to populate.
		fill(e.begin(), e.end(), 0);
		const size_t num_receptor_atoms = receptor_atoms.size();
		for (size_t l = 0; l < num_receptor_atoms; ++l)
		{
			const atom& a = atoms[receptor_atoms[l]];
			if (a.is_hydrogen()) continue;
			const double r2 = distance_sqr(probe_coords, a.coordinate);
			if (r2 <= scoring_function::cutoff_sqr)
			{
				const size_t t1 = a.xs;
				for (size_t i = 0; i < num_xs; ++i)
				{
					const size_t t2 = xs[i];
					const size_t type_pair_index = triangular_matrix_permissive_index(t1, t2);
					e[i] += sf.evaluate(type_pair_index, r2).e;
				}
			}
		}

		// Save accumulated free energies into grid maps.
		for (size_t i = 0; i < num_xs; ++i)
		{
			const size_t t = xs[i];
			grid_maps[t][grid_index(grid_idx)] = e[i];
		}
	}
}
