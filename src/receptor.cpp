#include <boost/filesystem/fstream.hpp>
#include "scoring_function.hpp"
#include "receptor.hpp"
#include "utility.hpp"

const float receptor::partition_granularity = 3.0f;
const float receptor::partition_granularity_inv = 1.0f / partition_granularity;

receptor::receptor(const path& p, const array<float, 3>& center, const array<float, 3>& size_, const float grid_granularity) : center(center), grid_granularity(grid_granularity), grid_granularity_inverse(1 / grid_granularity), grid_size(make_array(grid_granularity, grid_granularity, grid_granularity)), grid_size_inverse(make_array(grid_granularity_inverse, grid_granularity_inverse, grid_granularity_inverse)), grid_maps(scoring_function::n)
{
	// The loop may be unrolled by enabling compiler optimization.
	for (size_t i = 0; i < 3; ++i)
	{
		// Slightly expand the user-input size to the nearest multiple of granularity.
		num_grids[i] = static_cast<size_t>(ceil(size_[i] * grid_size_inverse[i]));
		size[i] = grid_size[i] * num_grids[i];
		num_probes[i] = num_grids[i] + 1;

		// Determine the two extreme corners.
		corner0[i] = center[i]  - size[i] * 0.5f;
		corner1[i] = corner0[i] + size[i];

		// Determine the number of partitions.
		num_partitions[i] = static_cast<size_t>(size[i] * partition_granularity_inv);
		partition_size[i] = size[i] / num_partitions[i];
		partition_size_inverse[i] = 1.0f / partition_size[i];
	}
	partitions.resize(num_partitions);

	// Parse the receptor line by line.
	atoms.reserve(2000); // A receptor typically consists of <= 2,000 atoms within bound.
	string residue = "XXXX"; // Current residue sequence located at 1-based [23, 26], used to track residue change, initialized to a dummy value.
	size_t residue_start; // The starting atom of the current residue.
	string line;
	for (boost::filesystem::ifstream ifs(p); getline(ifs, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "ATOM  " || record == "HETATM")
		{
			// If this line is the start of a new residue, mark the starting index within the atoms vector.
			if (line[25] != residue[3] || line[24] != residue[2] || line[23] != residue[1] || line[22] != residue[0])
			{
				residue[3] = line[25];
				residue[2] = line[24];
				residue[1] = line[23];
				residue[0] = line[22];
				residue_start = atoms.size();
			}

			// Parse the ATOM/HETATM line.
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
					if (b.is_hetero() && b.has_covalent_bond(a))
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
					if (!b.is_hetero() && b.has_covalent_bond(a))
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
					if (b.is_hetero() && b.has_covalent_bond(a))
					{
						a.dehydrophobicize();
						break;
					}
				}
			}
			if (project_distance_sqr(a.coord) < scoring_function::cutoff_sqr)
			{
				atoms.push_back(a);
			}
		}
		else if (record == "TER   ")
		{
			residue = "XXXX";
		}
	}

	// Allocate each nearby receptor atom to its corresponding partition.
	for (size_t z = 0; z < num_partitions[2]; ++z)
	for (size_t y = 0; y < num_partitions[1]; ++y)
	for (size_t x = 0; x < num_partitions[0]; ++x)
	{
		vector<size_t>& p = partitions(x, y, z);
		p.reserve(200);
		const array<size_t, 3> corner0_index = {x  , y  , z  };
		const array<size_t, 3> corner1_index = {x+1, y+1, z+1};
		const array<float, 3> corner0 = partition_corner0(corner0_index);
		const array<float, 3> corner1 = partition_corner0(corner1_index);
		for (size_t i = 0; i < atoms.size(); ++i)
		{
			if (project_distance_sqr(corner0, corner1, atoms[i].coord) < scoring_function::cutoff_sqr)
			{
				p.push_back(i);
			}
		}
	}
}

bool receptor::within(const array<float, 3>& coordinate) const
{
	for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
	{
		// Half-open-half-close box, i.e. [corner0, corner1)
		if (coordinate[i] < corner0[i] || corner1[i] <= coordinate[i])
			return false;
	}
	return true;
}

float receptor::project_distance_sqr(const array<float, 3>& corner0, const array<float, 3>& corner1, const array<float, 3>& coordinate) const
{
	// Calculate the projection point of the given coordinate onto the surface of the given box.
	array<float, 3> projection = coordinate; // The loop may be unrolled by enabling compiler optimization.
	for (size_t i = 0; i < 3; ++i)
	{
		if (projection[i] < corner0[i]) projection[i] = corner0[i];
		if (projection[i] > corner1[i]) projection[i] = corner1[i];
	}

	// Check if the distance between the projection and the given coordinate is within cutoff.
	return distance_sqr(projection, coordinate);
}

float receptor::project_distance_sqr(const array<float, 3>& coordinate) const
{
	return project_distance_sqr(corner0, corner1, coordinate);
}

array<float, 3> receptor::grid_corner0(const array<size_t, 3>& index) const
{
	return corner0 + (grid_size * index);
}

array<float, 3> receptor::partition_corner0(const array<size_t, 3>& index) const
{
	return corner0 + (partition_size * index);
}

array<size_t, 3> receptor::grid_index(const array<float, 3>& coordinate) const
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

array<size_t, 3> receptor::partition_index(const array<float, 3>& coordinate) const
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

int receptor::populate(const scoring_function& sf, const vector<size_t>& xs, const size_t z)
{
	const size_t n = xs.size();
	const size_t num_y_probes = num_probes[1];
	const size_t num_x_probes = num_probes[0];
	array<float, 3> probe_coord = corner0;
	probe_coord[2] += grid_granularity * z;
	size_t offset = num_x_probes * num_y_probes * z;
	vector<float> e(n);

	// For each probe atom coordinate of the given z dimension value.
	for (size_t y = 0; y < num_y_probes; ++y)
	{
		for (size_t x = 0; x < num_x_probes; ++x)
		{
			// Accumulate individual free energies for each atom types to populate.
			fill(e.begin(), e.end(), 0.0f);
			for (const auto p : partitions(partition_index(probe_coord)))
			{
				const atom& a = atoms[p];
				assert(!a.is_hydrogen());
				const float r2 = distance_sqr(probe_coord, a.coord);
				if (r2 < scoring_function::cutoff_sqr)
				{
					const size_t t1 = a.xs;
					for (size_t i = 0; i < n; ++i)
					{
						e[i] += sf.e[sf.o(mp(t1, xs[i]), r2)];
					}
				}
			}

			// Save accumulated free energies into grid maps.
			for (size_t i = 0; i < n; ++i)
			{
				grid_maps[xs[i]][offset] = e[i];
			}

			// Move to the next probe.
			++offset;
			probe_coord[0] += grid_granularity;
		}
		probe_coord[0] = corner0[0];
		probe_coord[1] += grid_granularity;
	}
	return 0;
}
