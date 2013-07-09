#include <iomanip>
#include <random>
#include "utility.hpp"
#include "ligand.hpp"

void solution::resize(const size_t nv, const size_t nf, const size_t na)
{
	x.resize(nv+1);
	a.resize(nf);
	q.resize(nf);
	c.resize(na);
	d.resize(na);
	f.resize(nf);
	t.resize(nf);
	g.resize(nv);
}

void frame::output(boost::filesystem::ofstream& ofs) const
{
	ofs << "BRANCH"    << setw(4) << rotorXsrn << setw(4) << rotorYsrn << '\n';
}

ligand::ligand(const path& p) : nt(0)
{
	// Initialize necessary variables for constructing a ligand.
	frames.reserve(30); // A ligand typically consists of <= 30 frames.
	frames.push_back(frame(0, 0, 1, 0, 0)); // ROOT is also treated as a frame. The parent and rotorX of ROOT frame are dummy.
	atoms.reserve(100); // A ligand typically consists of <= 100 heavy atoms.

	// Initialize helper variables for parsing.
	vector<vector<size_t>> bonds; ///< Covalent bonds.
	bonds.reserve(100); // A ligand typically consists of <= 100 heavy atoms.
	size_t current = 0; // Index of current frame, initialized to ROOT frame.
	frame* f = &frames.front(); // Pointer to the current frame.
	string line;

	// Parse the ligand line by line.
	for (boost::filesystem::ifstream ifs(p); getline(ifs, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "ATOM  " || record == "HETATM")
		{
			// Whenever an ATOM/HETATM line shows up, the current frame must be the last one.
			assert(current == frames.size() - 1);
			assert(f == &frames.back());

			// Parse the line.
			atom a(line);

			// Skip unsupported atom types.
			if (a.ad_unsupported()) continue;

			if (a.is_hydrogen()) // Current atom is a hydrogen.
			{
				for (size_t i = atoms.size(); i > f->beg;)
				{
					atom& b = atoms[--i];
					if (a.has_covalent_bond(b))
					{
						// For a polar hydrogen, the bonded hetero atom must be a hydrogen bond donor.
						if (a.is_polar_hydrogen())
						{
							assert(b.is_hetero());
							b.donorize();
						}
						else
						{
							assert(!b.is_hetero());
						}
						// Save the hydrogen.
						b.hydrogens.push_back(a);
						break;
					}
				}
			}
			else // Current atom is a heavy atom.
			{
				// Find bonds between the current atom and the other atoms of the same frame.
				assert(bonds.size() == atoms.size());
				bonds.push_back(vector<size_t>());
				bonds.back().reserve(4); // An atom typically consists of <= 4 bonds.
				for (size_t i = atoms.size(); i > f->beg;)
				{
					atom& b = atoms[--i];
					if (a.has_covalent_bond(b))
					{
						bonds[atoms.size()].push_back(i);
						bonds[i].push_back(atoms.size());

						// If carbon atom b is bonded to hetero atom a, b is no longer a hydrophobic atom.
						if (a.is_hetero() && !b.is_hetero())
						{
							b.dehydrophobicize();
						}
						// If carbon atom a is bonded to hetero atom b, a is no longer a hydrophobic atom.
						else if (!a.is_hetero() && b.is_hetero())
						{
							a.dehydrophobicize();
						}
					}
				}

				// Set rotorYidx if the serial number of current atom is rotorYsrn.
				if (a.serial == f->rotorYsrn)
				{
					f->rotorYidx = atoms.size();
					assert(f->rotorYidx == f->beg);
				}

				// Save the heavy atom.
				atoms.push_back(a);
			}
		}
		else if (record == "BRANCH")
		{
			// Parse "BRANCH   X   Y". X and Y are right-justified and 4 characters wide.
			const size_t rotorXsrn = stoul(line.substr( 6, 4));
			const size_t rotorYsrn = stoul(line.substr(10, 4));

			// Find the corresponding heavy atom with x as its atom serial number in the current frame.
			for (size_t i = f->beg; true; ++i)
			{
				if (atoms[i].serial == rotorXsrn)
				{
					// Insert a new frame whose parent is the current frame.
					frames.push_back(frame(current, rotorXsrn, rotorYsrn, i, atoms.size()));
					break;
				}
			}

			// The current frame has the newly inserted BRANCH frame as one of its branches.
			// It is unsafe to use f in place of frames[current] because frames could reserve a new memory block after calling push_back().
			frames[current].branches.push_back(frames.size() - 1);

			// Now the current frame is the newly inserted BRANCH frame.
			current = frames.size() - 1;

			// Update the pointer to the current frame.
			f = &frames[current];

			// The ending index of atoms of previous frame is the starting index of atoms of current frame.
			frames[current - 1].end = f->beg;
		}
		else if (record == "ENDBRA")
		{
			// A frame may be empty, e.g. "BRANCH   4   9" is immediately followed by "ENDBRANCH   4   9".
			// This emptiness is likely to be caused by invalid input structure, especially when all the atoms are located in the same plane.
			if (f->beg == atoms.size()) throw domain_error("Error parsing " + p.filename().string() + ": an empty BRANCH has been detected, indicating the input ligand structure is probably invalid.");

			// If the current frame consists of rotor Y and a few hydrogens only, e.g. -OH, -NH2 or -CH3,
			// the torsion of this frame will have no effect on scoring and is thus redundant.
			if (current + 1 == frames.size() && f->beg + 1 == atoms.size())
			{
				f->active = false;
			}
			else
			{
				++nt;
			}

			// Set up bonds between rotorX and rotorY.
			bonds[f->rotorYidx].push_back(f->rotorXidx);
			bonds[f->rotorXidx].push_back(f->rotorYidx);

			// Dehydrophobicize rotorX and rotorY if necessary.
			atom& rotorY = atoms[f->rotorYidx];
			atom& rotorX = atoms[f->rotorXidx];
			if (rotorY.is_hetero() && !rotorX.is_hetero()) rotorX.dehydrophobicize();
			if (rotorX.is_hetero() && !rotorY.is_hetero()) rotorY.dehydrophobicize();

			// Calculate yy and xy.
			const frame& p = frames[f->parent];
			f->yy = rotorY.coord - atoms[p.rotorYidx].coord;
			f->xy = normalize(rotorY.coord - rotorX.coord);

			// Now the parent of the following frame is the parent of current frame.
			current = f->parent;

			// Update the pointer to the current frame.
			f = &frames[current];
		}
	}
	assert(current == 0); // current should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.
	assert(f == &frames.front()); // The frame pointer should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.
	assert(frames.size() >= 1);
	assert(frames.size() - 1 >= nt);
	frames.back().end = atoms.size();
	nv = 6 + nt;

	// Update atoms[].coord relative to frame origin.
	for (const frame& f : frames)
	{
		const array<float, 3> origin = atoms[f.rotorYidx].coord;
		for (size_t i = f.beg; i < f.end; ++i)
		{
			atom& a = atoms[i];
			a.coord -= origin;
			for (atom& h : a.hydrogens)
			{
				h.coord -= origin;
			}
		}
	}

	// Find intra-ligand interacting pairs that are not 1-4.
	interacting_pairs.reserve(atoms.size() * atoms.size());
	vector<size_t> neighbors;
	neighbors.reserve(10); // An atom typically consists of <= 10 neighbors.
	for (size_t k1 = 0; k1 < frames.size(); ++k1)
	{
		const frame& f1 = frames[k1];
		for (size_t i = f1.beg; i < f1.end; ++i)
		{
			// Find neighbor atoms within 3 consecutive covalent bonds.
			for (const size_t b1 : bonds[i])
			{
				if (find(neighbors.cbegin(), neighbors.cend(), b1) == neighbors.cend())
				{
					neighbors.push_back(b1);
				}
				for (const size_t b2 : bonds[b1])
				{
					if (find(neighbors.cbegin(), neighbors.cend(), b2) == neighbors.cend())
					{
						neighbors.push_back(b2);
					}
					for (const size_t b3 : bonds[b2])
					{
						if (find(neighbors.cbegin(), neighbors.cend(), b3) == neighbors.cend())
						{
							neighbors.push_back(b3);
						}
					}
				}
			}

			// Determine if interacting pairs can be possibly formed.
			for (size_t k2 = k1 + 1; k2 < frames.size(); ++k2)
			{
				const frame& f2 = frames[k2];
				const frame& f3 = frames[f2.parent];
				for (size_t j = f2.beg; j < f2.end; ++j)
				{
					if (k1 == f2.parent && (i == f2.rotorXidx || j == f2.rotorYidx)) continue;
					if (k1 > 0 && f1.parent == f2.parent && i == f1.rotorYidx && j == f2.rotorYidx) continue;
					if (f2.parent > 0 && k1 == f3.parent && i == f3.rotorXidx && j == f2.rotorYidx) continue;
					if (find(neighbors.cbegin(), neighbors.cend(), j) != neighbors.cend()) continue;
					const size_t p_offset = scoring_function::nr * mp(atoms[i].xs, atoms[j].xs);
					interacting_pairs.push_back(interacting_pair(i, j, p_offset));
				}
			}

			// Clear the current neighbor set for the next atom.
			neighbors.clear();
		}
	}
}

bool ligand::evaluate(solution& s, const scoring_function& sf, const receptor& rec, const float e_upper_bound) const
{
	// Apply position and orientation to ROOT frame.
	s.c[0][0] = s.x[0];
	s.c[0][1] = s.x[1];
	s.c[0][2] = s.x[2];
	s.q[0][0] = s.x[3];
	s.q[0][1] = s.x[4];
	s.q[0][2] = s.x[5];
	s.q[0][3] = s.x[6];

	// Apply torsions to frames.
	for (size_t k = 0, t = 7; k < frames.size(); ++k)
	{
		const frame& f = frames[k];
		if (!f.active) continue;
		const float y_0 = s.c[f.rotorYidx][0];
		const float y_1 = s.c[f.rotorYidx][1];
		const float y_2 = s.c[f.rotorYidx][2];
		const float q_0 = s.q[k][0];
		const float q_1 = s.q[k][1];
		const float q_2 = s.q[k][2];
		const float q_3 = s.q[k][3];
		assert(fabs(q_0*q_0 + q_1*q_1 + q_2*q_2 + q_3*q_3 - 1.0f) < 1e-3f);
		const float q00 = q_0*q_0;
		const float q01 = q_0*q_1;
		const float q02 = q_0*q_2;
		const float q03 = q_0*q_3;
		const float q11 = q_1*q_1;
		const float q12 = q_1*q_2;
		const float q13 = q_1*q_3;
		const float q22 = q_2*q_2;
		const float q23 = q_2*q_3;
		const float q33 = q_3*q_3;
		const float m_0 = q00+q11-q22-q33;
		const float m_1 = 2*(-q03+q12);
		const float m_2 = 2*(q02+q13);
		const float m_3 = 2*(q03+q12);
		const float m_4 = q00-q11+q22-q33;
		const float m_5 = 2*(-q01+q23);
		const float m_6 = 2*(-q02+q13);
		const float m_7 = 2*(q01+q23);
		const float m_8 = q00-q11-q22+q33;
		for (size_t i = f.beg + 1; i < f.end; ++i)
		{
			const float o_0 = atoms[i].coord[0];
			const float o_1 = atoms[i].coord[1];
			const float o_2 = atoms[i].coord[2];
			s.c[i][0] = y_0 + m_0 * o_0 + m_1 * o_1 + m_2 * o_2;
			s.c[i][1] = y_1 + m_3 * o_0 + m_4 * o_1 + m_5 * o_2;
			s.c[i][2] = y_2 + m_6 * o_0 + m_7 * o_1 + m_8 * o_2;
		}
		for (const size_t i : f.branches)
		{
			const frame& b = frames[i];
			s.c[b.rotorYidx][0] = y_0 + m_0 * b.yy[0] + m_1 * b.yy[1] + m_2 * b.yy[2];
			s.c[b.rotorYidx][1] = y_1 + m_3 * b.yy[0] + m_4 * b.yy[1] + m_5 * b.yy[2];
			s.c[b.rotorYidx][2] = y_2 + m_6 * b.yy[0] + m_7 * b.yy[1] + m_8 * b.yy[2];

			// If the current BRANCH frame does not have an active torsion, skip it.
			if (!b.active)
			{
				assert(b.beg + 1 == b.end);
				assert(b.beg == b.rotorYidx);
				continue;
			}
			assert(normalized(b.xy));
			const float a_0 = m_0 * b.xy[0] + m_1 * b.xy[1] + m_2 * b.xy[2];
			const float a_1 = m_3 * b.xy[0] + m_4 * b.xy[1] + m_5 * b.xy[2];
			const float a_2 = m_6 * b.xy[0] + m_7 * b.xy[1] + m_8 * b.xy[2];
			assert(fabs(a_0*a_0 + a_1*a_1 + a_2*a_2 - 1.0f) < 1e-3f);
			s.a[i][0] = a_0;
			s.a[i][1] = a_1;
			s.a[i][2] = a_2;
			const float h = s.x[t++] * 0.5f;
			const float sinh = sin(h);
			const float r_0 = cos(h);
			const float r_1 = sinh * a_0;
			const float r_2 = sinh * a_1;
			const float r_3 = sinh * a_2;
			s.q[i][0] = r_0 * q_0 - r_1 * q_1 - r_2 * q_2 - r_3 * q_3;
			s.q[i][1] = r_0 * q_1 + r_1 * q_0 + r_2 * q_3 - r_3 * q_2;
			s.q[i][2] = r_0 * q_2 - r_1 * q_3 + r_2 * q_0 + r_3 * q_1;
			s.q[i][3] = r_0 * q_3 + r_1 * q_2 - r_2 * q_1 + r_3 * q_0;
			assert(normalized(s.q[i]));
		}
	}

	// Check steric clash between atoms of different frames except for (rotorX, rotorY) pair.
	//for (size_t k1 = frames.size() - 1; k1 > 0; --k1)
	//{
	//	const frame& f1 = frames[k1];
	//	for (size_t i1 = f1.beg; i1 < f1.end; ++i1)
	//	{
	//		for (size_t k2 = 0; k2 < k1; ++k2)
	//		{
	//			const frame& f2 = frames[k2];
	//			for (size_t i2 = f2.beg; i2 < f2.end; ++i2)
	//			{
	//				if ((distance_sqr(c[i1], c[i2]) < sqr(atoms[i1].covalent_radius() + atoms[i2].covalent_radius())) && (!((k2 == f1.parent) && (i1 == f1.rotorYidx) && (i2 == f1.rotorXidx))))
	//					return false;
	//			}
	//		}
	//	}
	//}

	float e = 0.0f;
	for (size_t i = 0; i < atoms.size(); ++i)
	{
		// Load coordinate from memory into registers.
		const float c_0 = s.c[i][0];
		const float c_1 = s.c[i][1];
		const float c_2 = s.c[i][2];

		// Deal with out-of-box case
		if (c_0 < rec.corner0[0] || rec.corner1[0] <= c_0 || c_1 < rec.corner0[1] || rec.corner1[1] <= c_1 || c_2 < rec.corner0[2] || rec.corner1[2] <= c_2)
		{
			e += 10;
			s.d[i][0] = 0;
			s.d[i][1] = 0;
			s.d[i][2] = 0;
			continue;
		}

		// Retrieve the grid map in need.
		const vector<float>& map = rec.maps[atoms[i].xs];
		assert(map.size());

		// Find the index of the current coordinates.
		const size_t index_0 = static_cast<size_t>((c_0 - rec.corner0[0]) * rec.granularity_inverse);
		const size_t index_1 = static_cast<size_t>((c_1 - rec.corner0[1]) * rec.granularity_inverse);
		const size_t index_2 = static_cast<size_t>((c_2 - rec.corner0[2]) * rec.granularity_inverse);
		assert(index_0 + 1 < rec.num_probes[0]);
		assert(index_1 + 1 < rec.num_probes[1]);
		assert(index_2 + 1 < rec.num_probes[2]);

		// Calculate the offsets to grid map and lookup the values.
		const size_t o000 = rec.num_probes[0] * (rec.num_probes[1] * index_2 + index_1) + index_0;
		const float e000 = map[o000];
		const float e100 = map[o000 + 1];
		const float e010 = map[o000 + rec.num_probes[0]];
		const float e001 = map[o000 + rec.num_probes[0] * rec.num_probes[1]];
		e += e000;
		s.d[i][0] = (e100 - e000) * rec.granularity_inverse;
		s.d[i][1] = (e010 - e000) * rec.granularity_inverse;
		s.d[i][2] = (e001 - e000) * rec.granularity_inverse;
	}

	// Calculate intra-ligand free energy.
	const size_t num_interacting_pairs = interacting_pairs.size();
	for (size_t i = 0; i < num_interacting_pairs; ++i)
	{
		const interacting_pair& p = interacting_pairs[i];
		const array<float, 3> r = s.c[p.i2] - s.c[p.i1];
		const float r2 = norm_sqr(r);
		if (r2 < scoring_function::cutoff_sqr)
		{
			const size_t o = p.p_offset + static_cast<size_t>(sf.ns * r2);
			e += sf.e[o];
			const array<float, 3> derivative = sf.d[o] * r;
			s.d[p.i1] -= derivative;
			s.d[p.i2] += derivative;
		}
	}

	// If the free energy is no better than the upper bound, refuse this conformation.
	if (e >= e_upper_bound) return false;

	// Store e from register into memory.
	s.e = e;

	// Calculate and aggregate the force and torque of BRANCH frames to their parent frame.
	fill(s.f.begin(), s.f.end(), zero3);
	fill(s.t.begin(), s.t.end(), zero3);
	for (size_t k = frames.size(), t = 6 + nt;;)
	{
		const frame& f = frames[--k];

		// Load variables from memory into registers.
		const float y_0 = s.c[f.rotorYidx][0];
		const float y_1 = s.c[f.rotorYidx][1];
		const float y_2 = s.c[f.rotorYidx][2];
		float f_0 = s.f[k][0];
		float f_1 = s.f[k][1];
		float f_2 = s.f[k][2];
		float t_0 = s.t[k][0];
		float t_1 = s.t[k][1];
		float t_2 = s.t[k][2];
		for (size_t i = f.beg; i < f.end; ++i)
		{
			const float d_0 = s.d[i][0];
			const float d_1 = s.d[i][1];
			const float d_2 = s.d[i][2];

			// The derivatives with respect to the position, orientation, and torsions
			// would be the negative total force acting on the ligand,
			// the negative total torque, and the negative torque projections, respectively,
			// where the projections refer to the torque applied to the branch moved by the torsion,
			// projected on its rotation axis.
			f_0 += d_0;
			f_1 += d_1;
			f_2 += d_2;
			const float yc_0 = s.c[i][0] - y_0;
			const float yc_1 = s.c[i][1] - y_1;
			const float yc_2 = s.c[i][2] - y_2;
			t_0 += yc_1 * d_2 - yc_2 * d_1;
			t_1 += yc_2 * d_0 - yc_0 * d_2;
			t_2 += yc_0 * d_1 - yc_1 * d_0;
		}

		// Save the aggregated force and torque of ROOT frame to g.
		if (k == 0)
		{
			s.g[0] = f_0;
			s.g[1] = f_1;
			s.g[2] = f_2;
			s.g[3] = t_0;
			s.g[4] = t_1;
			s.g[5] = t_2;
			return true;
		}

		// Store variables from registers into memory.
		s.f[k][0] = f_0;
		s.f[k][1] = f_1;
		s.f[k][2] = f_2;
		s.t[k][0] = t_0;
		s.t[k][1] = t_1;
		s.t[k][2] = t_2;

		// Aggregate the force and torque of current frame to its parent frame.
		s.f[f.parent][0] += f_0;
		s.f[f.parent][1] += f_1;
		s.f[f.parent][2] += f_2;
		const size_t p = frames[f.parent].rotorYidx;
		const float yy_0 = y_0 - s.c[p][0];
		const float yy_1 = y_1 - s.c[p][1];
		const float yy_2 = y_2 - s.c[p][2];
		s.t[f.parent][0] += t_0 + yy_1 * f_2 - yy_2 * f_1;
		s.t[f.parent][1] += t_1 + yy_2 * f_0 - yy_0 * f_2;
		s.t[f.parent][2] += t_2 + yy_0 * f_1 - yy_1 * f_0;

		// If the current BRANCH frame does not have an active torsion, skip it.
		if (!f.active) continue;

		// Save the aggregated torque of BRANCH frames to g.
		s.g[--t] = t_0 * s.a[k][0] + t_1 * s.a[k][1] + t_2 * s.a[k][2]; // dot product
	}
}

int ligand::bfgs(solution& s0, const scoring_function& sf, const receptor& rec, const size_t seed, const size_t num_generations) const
{
	const size_t num_alphas = 5; // Number of alpha values for determining step size in BFGS
	const float e_upper_bound = 40.0f * atoms.size(); // A conformation will be droped if its free energy is not better than e_upper_bound.
	solution s1, s2;
	s0.resize(nv, frames.size(), atoms.size());
	s1.resize(nv, frames.size(), atoms.size());
	s2.resize(nv, frames.size(), atoms.size());
	vector<float> p(nv), y(nv), mhy(nv);
	vector<float> h(nv*(nv+1)>>1); // Symmetric triangular Hessian matrix.
	float alpha, pg1, pg2, yhy, yp, ryp, pco, qw, qx, qy, qz, qnorm_inv;
	size_t g, i, j;
	mt19937_64 rng(seed);
	uniform_real_distribution<float> uniform_11(-1.0f, 1.0f);

	// Randomize conformation x0.
	s0.x[0] = rec.center[0] + uniform_11(rng) * rec.size[0];
	s0.x[1] = rec.center[1] + uniform_11(rng) * rec.size[1];
	s0.x[2] = rec.center[2] + uniform_11(rng) * rec.size[2];
	qw = uniform_11(rng);
	qx = uniform_11(rng);
	qy = uniform_11(rng);
	qz = uniform_11(rng);
	qnorm_inv = 1.0f / sqrt(qw*qw + qx*qx + qy*qy + qz*qz);
	s0.x[3] = qw * qnorm_inv;
	s0.x[4] = qx * qnorm_inv;
	s0.x[5] = qy * qnorm_inv;
	s0.x[6] = qz * qnorm_inv;
	for (i = 0; i < nt; ++i)
	{
		s0.x[7 + i] = uniform_11(rng);
	}
	evaluate(s0, sf, rec, e_upper_bound);

	for (g = 0; g < num_generations; ++g)
	{
		// Make a copy, so the previous conformation is retained.
		s1.x = s0.x;
		s1.x[0] += uniform_11(rng);
		s1.x[1] += uniform_11(rng);
		s1.x[2] += uniform_11(rng);
		evaluate(s1, sf, rec, e_upper_bound);

		// Initialize the inverse Hessian matrix to identity matrix.
		// An easier option that works fine in practice is to use a scalar multiple of the identity matrix,
		// where the scaling factor is chosen to be in the range of the eigenvalues of the true Hessian.
		// See N&R for a recipe to find this initializer.
		fill(h.begin(), h.end(), 0.0f);
		for (i = 0; i < nv; ++i)
			h[mr(i, i)] = 1.0f;

		// Given the mutated conformation c1, use BFGS to find a local minimum.
		// The conformation of the local minimum is saved to c2, and its derivative is saved to g2.
		// http://en.wikipedia.org/wiki/BFGS_method
		// http://en.wikipedia.org/wiki/Quasi-Newton_method
		// The loop breaks when an appropriate alpha cannot be found.
		while (true)
		{
			// Calculate p = -h*g, where p is for descent direction, h for Hessian, and g for gradient.
			for (i = 0; i < nv; ++i)
			{
				float sum = 0.0f;
				for (j = 0; j < nv; ++j)
					sum += h[mp(i, j)] * s1.g[j];
				p[i] = -sum;
			}

			// Calculate pg = p*g = -h*g^2 < 0
			pg1 = 0;
			for (i = 0; i < nv; ++i)
				pg1 += p[i] * s1.g[i];

			// Perform a line search to find an appropriate alpha.
			// Try different alpha values for num_alphas times.
			// alpha starts with 1, and shrinks to alpha_factor of itself iteration by iteration.
			alpha = 1.0;
			for (j = 0; j < num_alphas; ++j)
			{
				// Calculate c2 = c1 + ap.
				s2.x[0] = s1.x[0] + alpha * p[0];
				s2.x[1] = s1.x[1] + alpha * p[1];
				s2.x[2] = s1.x[2] + alpha * p[2];
				const float po0 = p[3];
				const float po1 = p[4];
				const float po2 = p[5];
				const float ponrm = sqrt(po0*po0 + po1*po1 + po2*po2);
				const float h = 0.5f * alpha * ponrm;
				const float u = sin(h) / ponrm;
				const float pq0 = cos(h);
				const float pq1 = u*po0;
				const float pq2 = u*po1;
				const float pq3 = u*po2;
				assert(fabs(pq0*pq0 + pq1*pq1 + pq2*pq2 + pq3*pq3 - 1.0f) < 1e-3f);
				const float x1q0 = s1.x[3];
				const float x1q1 = s1.x[4];
				const float x1q2 = s1.x[5];
				const float x1q3 = s1.x[6];
				assert(fabs(x1q0*x1q0 + x1q1*x1q1 + x1q2*x1q2 + x1q3*x1q3 - 1.0f) < 1e-3f);
				const float x2q0 = pq0 * x1q0 - pq1 * x1q1 - pq2 * x1q2 - pq3 * x1q3;
				const float x2q1 = pq0 * x1q1 + pq1 * x1q0 + pq2 * x1q3 - pq3 * x1q2;
				const float x2q2 = pq0 * x1q2 - pq1 * x1q3 + pq2 * x1q0 + pq3 * x1q1;
				const float x2q3 = pq0 * x1q3 + pq1 * x1q2 - pq2 * x1q1 + pq3 * x1q0;
				assert(fabs(x2q0*x2q0 + x2q1*x2q1 + x2q2*x2q2 + x2q3*x2q3 - 1.0f) < 1e-3f);
				s2.x[3] = x2q0;
				s2.x[4] = x2q1;
				s2.x[5] = x2q2;
				s2.x[6] = x2q3;
				for (i = 0; i < nt; ++i)
				{
					s2.x[7 + i] = s1.x[7 + i] + alpha * p[6 + i];
				}

				// Evaluate c2, subject to Wolfe conditions http://en.wikipedia.org/wiki/Wolfe_conditions
				// 1) Armijo rule ensures that the step length alpha decreases f sufficiently.
				// 2) The curvature condition ensures that the slope has been reduced sufficiently.
				if (evaluate(s2, sf, rec, s1.e + 0.0001f * alpha * pg1))
				{
					pg2 = 0;
					for (i = 0; i < nv; ++i)
						pg2 += p[i] * s2.g[i];
					if (pg2 >= 0.9f * pg1)
						break; // An appropriate alpha is found.
				}

				alpha *= 0.1f;
			}

			// If an appropriate alpha cannot be found, exit the BFGS loop.
			if (j == num_alphas) break;

			// Update Hessian matrix h.
			for (i = 0; i < nv; ++i) // Calculate y = g2 - g1.
				y[i] = s2.g[i] - s1.g[i];
			for (i = 0; i < nv; ++i) // Calculate mhy = -h * y.
			{
				float sum = 0.0f;
				for (j = 0; j < nv; ++j)
					sum += h[mp(i, j)] * y[j];
				mhy[i] = -sum;
			}
			yhy = 0;
			for (i = 0; i < nv; ++i) // Calculate yhy = -y * mhy = -y * (-hy).
				yhy -= y[i] * mhy[i];
			yp = 0;
			for (i = 0; i < nv; ++i) // Calculate yp = y * p.
				yp += y[i] * p[i];
			ryp = 1 / yp;
			pco = ryp * (ryp * yhy + alpha);
			for (i = 0; i < nv; ++i)
			for (j = i; j < nv; ++j) // includes i
			{
				h[mr(i, j)] += ryp * (mhy[i] * p[j] + mhy[j] * p[i]) + pco * p[i] * p[j];
			}

			// Move to the next iteration.
			s1.x = s2.x;
			s1.e = s2.e;
			s1.g = s2.g;
		}

		// Accept c1 according to Metropolis criteria.
		if (s1.e < s0.e)
		{
			s0.x = s1.x;
			s0.e = s1.e;
		}
	}
	return 0;
}

void ligand::save(const path& output_ligand_path, const ptr_vector<solution>& solutions, const vector<size_t>& representatives) const
{
	assert(representatives.size());
	assert(representatives.size() <= solutions.size());
	boost::filesystem::ofstream ofs(output_ligand_path);
	ofs.setf(ios::fixed, ios::floatfield);
	ofs << setprecision(3);
	for (size_t k = 0; k < representatives.size(); ++k)
	{
		const solution& s = solutions[representatives[k]];
		ofs << "MODEL     " << setw(4) << (k + 1) << '\n'
			<< "REMARK     pKd:" << setw(7) << s.e << '\n';

		// Dump the ROOT frame.
		ofs << "ROOT\n";
		{
			const frame& f = frames.front();
			const array<float, 9> m = qtn4_to_mat3(s.q.front());
			for (size_t i = f.beg; i < f.end; ++i)
			{
				const atom& a = atoms[i];
				a.output(ofs, s.c[i]);
				if (a.hydrogens.empty()) continue;
				for (const atom& h : a.hydrogens)
				{
					h.output(ofs, s.c[f.rotorYidx] + m * h.coord);
				}
			}
		}
		ofs << "ENDROOT\n";

		// Dump the BRANCH frames.
		vector<bool> dumped(frames.size()); // dump_branches[0] is dummy. The ROOT frame has been dumped.
		vector<size_t> stack; // Stack to track the depth-first traversal sequence of frames in order to avoid recursion.
		stack.reserve(frames.size() - 1); // The ROOT frame is excluded.
		{
			const frame& f = frames.front();
			for (auto i = f.branches.rbegin(); i < f.branches.rend(); ++i)
			{
				stack.push_back(*i);
			}
		}
		while (!stack.empty())
		{
			const size_t fn = stack.back();
			const frame& f = frames[fn];
			if (dumped[fn]) // This BRANCH frame has been dumped.
			{
				ofs << "END";
				f.output(ofs);
				stack.pop_back();
			}
			else // This BRANCH frame has not been dumped.
			{
				f.output(ofs);
				const array<float, 9> m = qtn4_to_mat3(f.active ? s.q[fn] : s.q[f.parent]);
				for (size_t i = f.beg; i < f.end; ++i)
				{
					const atom& a = atoms[i];
					a.output(ofs, s.c[i]);
					if (a.hydrogens.empty()) continue;
					for (const atom& h : a.hydrogens)
					{
						h.output(ofs, s.c[f.rotorYidx] + m * h.coord);
					}
				}
				dumped[fn] = true;
				for (auto i = f.branches.rbegin(); i < f.branches.rend(); ++i)
				{
					stack.push_back(*i);
				}
			}
		}
		ofs << "TORSDOF " << frames.size() - 1 << '\n'
			<< "ENDMDL" << '\n';
	}
}
