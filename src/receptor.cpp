#include <cmath>
#include <cassert>
#include <boost/filesystem/fstream.hpp>
#include "scoring_function.hpp"
#include "array.hpp"
#include "receptor.hpp"

const double receptor::Default_Partition_Granularity = static_cast<double>(3);
const double receptor::Default_Partition_Granularity_Inverse = 1 / Default_Partition_Granularity;

receptor::receptor(const path& p, const array<double, 3>& center, const array<double, 3>& span_, const double grid_granularity) : center(center), grid_granularity(grid_granularity), grid_granularity_inverse(1 / grid_granularity), grid_size({ grid_granularity, grid_granularity, grid_granularity }), grid_size_inverse({ grid_granularity_inverse, grid_granularity_inverse, grid_granularity_inverse }), grid_maps(scoring_function::n)
{
	// The loop may be unrolled by enabling compiler optimization.
	for (size_t i = 0; i < 3; ++i)
	{
		// Slightly expand the user-input span to the nearest multiple of granularity.
		num_grids[i] = static_cast<size_t>(ceil(span_[i] * grid_size_inverse[i]));
		assert(num_grids[i] > 0);
		span[i] = grid_size[i] * num_grids[i];
		num_probes[i] = num_grids[i] + 1;

		// Determine the two extreme corners.
		corner1[i] = center[i] - span[i] * static_cast<double>(0.5);
		corner2[i] = corner1[i] + span[i];
		assert(corner1[i] < corner2[i]);

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
			}
			else if (a.is_hetero()) // It is a hetero atom.
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
			else // It is a carbon atom.
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
			atoms.push_back(a);
		}
		else if (record == "TER   ")
		{
			residue = "XXXX";
		}
	}

	// Find all the heavy receptor atoms that are within 8A of the box.
	vector<size_t> receptor_atoms_within_cutoff;
	receptor_atoms_within_cutoff.reserve(atoms.size());
	const size_t num_rec_atoms = atoms.size();
	for (size_t i = 0; i < num_rec_atoms; ++i)
	{
		const atom& a = atoms[i];
		if (project_distance_sqr(a.coordinate) < scoring_function::Cutoff_Sqr)
		{
			receptor_atoms_within_cutoff.push_back(i);
		}
	}
	const size_t num_receptor_atoms_within_cutoff = receptor_atoms_within_cutoff.size();

	// Allocate each nearby receptor atom to its corresponding partition.
	for (size_t x = 0; x < num_partitions[0]; ++x)
	for (size_t y = 0; y < num_partitions[1]; ++y)
	for (size_t z = 0; z < num_partitions[2]; ++z)
	{
		vector<size_t>& par = partitions[partition_index(array<size_t, 3>{x, y, z})];
		par.reserve(num_receptor_atoms_within_cutoff);
		const array<size_t, 3> index1 = {{ x,     y,     z     }};
		const array<size_t, 3> index2 = {{ x + 1, y + 1, z + 1 }};
		const array<double, 3> corner1 = partition_corner1(index1);
		const array<double, 3> corner2 = partition_corner1(index2);
		for (size_t l = 0; l < num_receptor_atoms_within_cutoff; ++l)
		{
			const size_t i = receptor_atoms_within_cutoff[l];
			const atom& a = atoms[i];
			const double proj_dist_sqr = project_distance_sqr(corner1, corner2, a.coordinate);
			if (proj_dist_sqr < scoring_function::Cutoff_Sqr)
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
		// Half-open-half-close box, i.e. [corner1, corner2)
		if (coordinate[i] < corner1[i] || corner2[i] <= coordinate[i])
			return false;
	}
	return true;
}

double receptor::project_distance_sqr(const array<double, 3>& corner1, const array<double, 3>& corner2, const array<double, 3>& coordinate) const
{
	// Calculate the projection point of the given coordinate onto the surface of the given box.
	array<double, 3> projection = coordinate; // The loop may be unrolled by enabling compiler optimization.
	for (size_t i = 0; i < 3; ++i)
	{
		if (projection[i] < corner1[i]) projection[i] = corner1[i];
		if (projection[i] > corner2[i]) projection[i] = corner2[i];
	}

	// Check if the distance between the projection and the given coordinate is within cutoff.
	return distance_sqr(projection, coordinate);
}

double receptor::project_distance_sqr(const array<double, 3>& coordinate) const
{
	return project_distance_sqr(corner1, corner2, coordinate);
}

array<double, 3> receptor::grid_corner1(const array<size_t, 3>& index) const
{
	return corner1 + (grid_size * index);
}

array<double, 3> receptor::partition_corner1(const array<size_t, 3>& index) const
{
	return corner1 + (partition_size * index);
}

array<size_t, 3> receptor::grid_index(const array<double, 3>& coordinate) const
{
	array<size_t, 3> index;
	for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
	{
		index[i] = static_cast<size_t>((coordinate[i] - corner1[i]) * grid_size_inverse[i]);
		// Boundary checking is not necessary because the given coordinate is a ligand atom,
		// which has been restricted within the half-open-half-close box [corner1, corner2).
		//if (index[i] == num_grids[i]) index[i] = num_grids[i] - 1;
	}
	return index;
}

array<size_t, 3> receptor::partition_index(const array<double, 3>& coordinate) const
{
	array<size_t, 3> index;
	for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
	{
		index[i] = static_cast<size_t>((coordinate[i] - corner1[i]) * partition_size_inverse[i]);
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
