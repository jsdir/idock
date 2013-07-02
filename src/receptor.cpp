#include <boost/filesystem/fstream.hpp>
#include "scoring_function.hpp"
#include "receptor.hpp"
#include "utility.hpp"

receptor::receptor(const path& p, const array<float, 3>& center, const array<float, 3>& size, const float granularity) : center(center), size(size), corner0(center - 0.5f * size), corner1(corner0 + size), granularity(granularity), granularity_inverse(1.0f / granularity), num_probes_product(1), p_offset(scoring_function::n), maps(scoring_function::n)
{
	// The loop may be unrolled by enabling compiler optimization.
	for (size_t i = 0; i < 3; ++i)
	{
		// Reserve one more probe to calculate the derivative.
		num_probes[i] = static_cast<size_t>(size[i] * granularity_inverse) + 1;
		num_probes_product *= num_probes[i];
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

			// Save the atom if and only if its distance to its projection point on the box is within cutoff.
			float r2 = 0.0f;
			for (size_t i = 0; i < 3; ++i)
			{
				if (a.coord[i] < corner0[i])
				{
					const float d = a.coord[i] - corner0[i];
					r2 += d * d;
				}
				else if (a.coord[i] > corner1[i])
				{
					const float d = a.coord[i] - corner1[i];
					r2 += d * d;
				}
			}
			if (r2 < scoring_function::cutoff_sqr)
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

int receptor::populate(const scoring_function& sf, const vector<size_t>& xs, const size_t z)
{
	const size_t n = xs.size();
	const float z_coord = corner0[2] + granularity * z;
	const size_t z_offset = num_probes[0] * num_probes[1] * z;

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
		const size_t y_beg = y_lb > corner0[1] ? (y_lb < corner1[1] ? static_cast<size_t>((y_lb - corner0[1]) * granularity_inverse)     : num_probes[1]) : 0;
		const size_t y_end = y_ub > corner0[1] ? (y_ub < corner1[1] ? static_cast<size_t>((y_ub - corner0[1]) * granularity_inverse) + 1 : num_probes[1]) : 0;
		const vector<size_t>& p = p_offset[a.xs];
		size_t zy_offset = z_offset + num_probes[0] * y_beg;
		float dy = corner0[1] + granularity * y_beg - a.coord[1];
		for (size_t y = y_beg; y < y_end; ++y, zy_offset += num_probes[0], dy += granularity)
		{
			const float dy_sqr = dy * dy;
			const float dx_sqr_ub = dydx_sqr_ub - dy_sqr;
			if (dx_sqr_ub <= 0) continue;
			const float dx_ub = sqrt(dx_sqr_ub);
			const float x_lb = a.coord[0] - dx_ub;
			const float x_ub = a.coord[0] + dx_ub;
			const size_t x_beg = x_lb > corner0[0] ? (x_lb < corner1[0] ? static_cast<size_t>((x_lb - corner0[0]) * granularity_inverse)     : num_probes[0]) : 0;
			const size_t x_end = x_ub > corner0[0] ? (x_ub < corner1[0] ? static_cast<size_t>((x_ub - corner0[0]) * granularity_inverse) + 1 : num_probes[0]) : 0;
			const float dzdy_sqr = dz_sqr + dy_sqr;
			size_t zyx_offset = zy_offset + x_beg;
			float dx = corner0[0] + granularity * x_beg - a.coord[0];
			for (size_t x = x_beg; x < x_end; ++x, ++zyx_offset, dx += granularity)
			{
				const float dx_sqr = dx * dx;
				const float r2 = dzdy_sqr + dx_sqr;
				if (r2 >= scoring_function::cutoff_sqr) continue;
				const size_t r_offset = static_cast<size_t>(sf.ns * r2);
				for (size_t i = 0; i < n; ++i)
				{
					maps[xs[i]][zyx_offset] += sf.e[p[i] + r_offset];
				}
			}
		}
	}
	return 0;
}
