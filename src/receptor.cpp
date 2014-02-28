#include <boost/filesystem/fstream.hpp>
#include "scoring_function.hpp"
#include "receptor.hpp"

receptor::receptor(const path& p, const box& b) : partitions(b.num_partitions[0] * b.num_partitions[1] * b.num_partitions[2]), grid_maps(scoring_function::n)
{
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
		if (b.project_distance_sqr(a.coordinate) < scoring_function::Cutoff_Sqr)
		{
			receptor_atoms_within_cutoff.push_back(i);
		}
	}
	const size_t num_receptor_atoms_within_cutoff = receptor_atoms_within_cutoff.size();

	// Allocate each nearby receptor atom to its corresponding partition.
	for (size_t x = 0; x < b.num_partitions[0]; ++x)
	for (size_t y = 0; y < b.num_partitions[1]; ++y)
	for (size_t z = 0; z < b.num_partitions[2]; ++z)
	{
		vector<size_t>& par = partitions[b.partition_index(array<size_t, 3>{x, y, z})];
		par.reserve(num_receptor_atoms_within_cutoff);
		const array<size_t, 3> index1 = {{ x,     y,     z     }};
		const array<size_t, 3> index2 = {{ x + 1, y + 1, z + 1 }};
		const array<double, 3> corner1 = b.partition_corner1(index1);
		const array<double, 3> corner2 = b.partition_corner1(index2);
		for (size_t l = 0; l < num_receptor_atoms_within_cutoff; ++l)
		{
			const size_t i = receptor_atoms_within_cutoff[l];
			const atom& a = atoms[i];
			const double proj_dist_sqr = b.project_distance_sqr(corner1, corner2, a.coordinate);
			if (proj_dist_sqr < scoring_function::Cutoff_Sqr)
			{
				par.push_back(i);
			}
		}
	}
}
