#include <boost/filesystem/fstream.hpp>
#include "scoring_function.hpp"
#include "receptor.hpp"
#include "utility.hpp"

receptor::receptor(const path& p, const array<float, 3>& center, const array<float, 3>& size_, const float granularity) : center(center), granularity(granularity), granularity_inverse(1.0f / granularity), grid_maps(scoring_function::n)
{
	// The loop may be unrolled by enabling compiler optimization.
	for (size_t i = 0; i < 3; ++i)
	{
		// Slightly expand the user-input size to the nearest multiple of granularity.
		num_grids[i] = static_cast<size_t>(ceil(size_[i] * granularity_inverse));
		size[i] = granularity * num_grids[i];
		num_probes[i] = num_grids[i] + 1;

		// Determine the two extreme corners.
		corner0[i] = center[i]  - size[i] * 0.5f;
		corner1[i] = corner0[i] + size[i];
	}

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

array<size_t, 3> receptor::grid_index(const array<float, 3>& coordinate) const
{
	array<size_t, 3> index;
	for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
	{
		index[i] = static_cast<size_t>((coordinate[i] - corner0[i]) * granularity_inverse);
		// Boundary checking is not necessary because the given coordinate is a ligand atom,
		// which has been restricted within the half-open-half-close box [corner0, corner1).
		//if (index[i] == num_grids[i]) index[i] = num_grids[i] - 1;
	}
	return index;
}

int receptor::populate(const scoring_function& sf, const vector<size_t>& xs, const size_t z)
{
	const size_t n = xs.size();
	const size_t num_y_probes = num_probes[1];
	const size_t num_x_probes = num_probes[0];
	const float z_coord = corner0[2] + granularity * z;
	const size_t z_offset = num_x_probes * num_y_probes * z;

	for (const auto& a : atoms)
	{
		assert(!a.is_hydrogen());
		const float dz = z_coord - a.coord[2];
		const float dz_sqr = dz * dz;
		const float dydx_sqr_ub = scoring_function::cutoff_sqr - dz_sqr;
		if (dydx_sqr_ub <= 0) continue;
		const float dydx_ub = sqrt(dydx_sqr_ub);
		const float y_lb = a.coord[1] - dydx_ub;
		const float y_ub = a.coord[1] + dydx_ub;
		const int y_begin = y_lb > corner0[1] ? static_cast<int>((y_lb - corner0[1]) * granularity_inverse) : 0;
		const int y_end = y_ub < corner1[1] ? static_cast<int>((y_ub - corner0[1]) * granularity_inverse) : num_y_probes;
		const size_t t1 = a.xs;
		size_t zy_offset = z_offset + num_x_probes * y_begin;
		float dy = corner0[1] + granularity * y_begin - a.coord[1];
		for (int y = y_begin; y < y_end; ++y, zy_offset += num_x_probes, dy += granularity)
		{
			const float dy_sqr = dy * dy;
			const float dx_sqr_ub = dydx_sqr_ub - dy_sqr;
			if (dx_sqr_ub <= 0) continue;
			const float dx_ub = sqrt(dx_sqr_ub);
			const float x_lb = a.coord[0] - dx_ub;
			const float x_ub = a.coord[0] + dx_ub;
			const int x_begin = x_lb > corner0[0] ? static_cast<int>((x_lb - corner0[0]) * granularity_inverse) : 0;
			const int x_end = x_ub < corner1[0] ? static_cast<int>((x_ub - corner0[0]) * granularity_inverse) : num_x_probes;
			const float dzdy_sqr = dz_sqr + dy_sqr;
			size_t zyx_offset = zy_offset + x_begin;
			float dx = corner0[0] + granularity * x_begin - a.coord[0];
			for (int x = x_begin; x < x_end; ++x, ++zyx_offset, dx += granularity)
			{
				const float dx_sqr = dx * dx;
				const float r2 = dzdy_sqr + dx_sqr;
				if (r2 >= scoring_function::cutoff_sqr) continue;
				for (size_t i = 0; i < n; ++i)
				{
					const size_t t2 = xs[i];
					grid_maps[t2][zyx_offset] += sf.e[sf.o(mp(t1, t2), r2)];
				}
			}
		}
	}
	return 0;
}
