#include <iomanip>
#include <random>
#include "utility.hpp"
#include "ligand.hpp"

void solution::resize(const size_t nv, const size_t nf, const size_t na)
{
	x.resize(nv+1);
	a.resize(3 * nf);
	q.resize(nf);
	c.resize(na);
	d.resize(na);
	f.resize(3 * nf);
	t.resize(3 * nf);
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
	float e, y0, y1, y2, q0, q1, q2, q3, q00, q01, q02, q03, q11, q12, q13, q22, q23, q33, m0, m1, m2, m3, m4, m5, m6, m7, m8, v0, v1, v2, c0, c1, c2, e000, e100, e010, e001, a0, a1, a2, h, sinh, r0, r1, r2, r3, vs, dor, f0, f1, f2, t0, t1, t2, d0, d1, d2;
	size_t k, t, i, i0, i1, i2, o0, o1, o2;

	// Apply position, orientation and torsions.
	s.c[0][0] = s.x[0];
	s.c[0][1] = s.x[1];
	s.c[0][2] = s.x[2];
	s.q[0][0] = s.x[3];
	s.q[0][1] = s.x[4];
	s.q[0][2] = s.x[5];
	s.q[0][3] = s.x[6];
	e = 0.0f;
	for (k = 0, t = 7; k < frames.size(); ++k)
	{
		const frame& f = frames[k];
		y0 = s.c[f.rotorYidx][0];
		y1 = s.c[f.rotorYidx][1];
		y2 = s.c[f.rotorYidx][2];
		if (f.active)
		{
			// Translate orientation from quaternion into 3x3 matrix.
			q0 = s.q[k][0];
			q1 = s.q[k][1];
			q2 = s.q[k][2];
			q3 = s.q[k][3];
			assert(fabs(q0*q0 + q1*q1 + q2*q2 + q3*q3 - 1.0f) < 1e-3f);
			q00 = q0 * q0;
			q01 = q0 * q1;
			q02 = q0 * q2;
			q03 = q0 * q3;
			q11 = q1 * q1;
			q12 = q1 * q2;
			q13 = q1 * q3;
			q22 = q2 * q2;
			q23 = q2 * q3;
			q33 = q3 * q3;
			m0 = q00 + q11 - q22 - q33;
			m1 = 2 * (q12 - q03);
			m2 = 2 * (q02 + q13);
			m3 = 2 * (q03 + q12);
			m4 = q00 - q11 + q22 - q33;
			m5 = 2 * (q23 - q01);
			m6 = 2 * (q13 - q02);
			m7 = 2 * (q01 + q23);
			m8 = q00 - q11 - q22 + q33;
		}
		for (i = f.beg; i < f.end; ++i)
		{
			const atom& a = atoms[i];
			if (i == f.beg)
			{
				c0 = y0;
				c1 = y1;
				c2 = y2;
			}
			else
			{
				// Calculate coordinate from transformation matrix and offset.
				v0 = a.coord[0];
				v1 = a.coord[1];
				v2 = a.coord[2];
				c0 = y0 + m0 * v0 + m1 * v1 + m2 * v2;
				c1 = y1 + m3 * v0 + m4 * v1 + m5 * v2;
				c2 = y2 + m6 * v0 + m7 * v1 + m8 * v2;

				// Store coordinate from registers into memory.
				s.c[i][0] = c0;
				s.c[i][1] = c1;
				s.c[i][2] = c2;
			}
			// Penalize out-of-box case.
			if (c0 < rec.corner0[0] || rec.corner1[0] <= c0 || c1 < rec.corner0[1] || rec.corner1[1] <= c1 || c2 < rec.corner0[2] || rec.corner1[2] <= c2)
			{
				e += 10;
				s.d[i][0] = 0;
				s.d[i][1] = 0;
				s.d[i][2] = 0;
				continue;
			}

			// Find the index of the current coordinates.
			i0 = static_cast<size_t>((c0 - rec.corner0[0]) * rec.granularity_inverse);
			i1 = static_cast<size_t>((c1 - rec.corner0[1]) * rec.granularity_inverse);
			i2 = static_cast<size_t>((c2 - rec.corner0[2]) * rec.granularity_inverse);
			assert(i0 + 1 < rec.num_probes[0]);
			assert(i1 + 1 < rec.num_probes[1]);
			assert(i2 + 1 < rec.num_probes[2]);
			o0 = rec.num_probes[0] * (rec.num_probes[1] * i2 + i1) + i0;

			// Retrieve the grid map and lookup the values.
			const vector<float>& map = rec.maps[a.xs];
			assert(map.size());
			e000 = map[o0];
			e100 = map[o0 + 1];
			e010 = map[o0 + rec.num_probes[0]];
			e001 = map[o0 + rec.num_probes[0] * rec.num_probes[1]];
			e += e000;
			s.d[i][0] = (e100 - e000) * rec.granularity_inverse;
			s.d[i][1] = (e010 - e000) * rec.granularity_inverse;
			s.d[i][2] = (e001 - e000) * rec.granularity_inverse;
		}
		for (const size_t i : f.branches)
		{
			const frame& b = frames[i];
			s.c[b.rotorYidx][0] = y0 + m0 * b.yy[0] + m1 * b.yy[1] + m2 * b.yy[2];
			s.c[b.rotorYidx][1] = y1 + m3 * b.yy[0] + m4 * b.yy[1] + m5 * b.yy[2];
			s.c[b.rotorYidx][2] = y2 + m6 * b.yy[0] + m7 * b.yy[1] + m8 * b.yy[2];

			// Skip inactive BRANCH frames.
			if (!b.active) continue;

			// Update s.a of BRANCH frames.
			a0 = m0 * b.xy[0] + m1 * b.xy[1] + m2 * b.xy[2];
			a1 = m3 * b.xy[0] + m4 * b.xy[1] + m5 * b.xy[2];
			a2 = m6 * b.xy[0] + m7 * b.xy[1] + m8 * b.xy[2];
			assert(fabs(a0*a0 + a1*a1 + a2*a2 - 1.0f) < 1e-3f);
			o0 = 3 * i;
			o1 = o0 + 1;
			o2 = o1 + 1;
			s.a[o0] = a0;
			s.a[o1] = a1;
			s.a[o2] = a2;

			// Update s.q of BRANCH frames.
			h = s.x[t++] * 0.5f;
			sinh = sin(h);
			r0 = cos(h);
			r1 = sinh * a0;
			r2 = sinh * a1;
			r3 = sinh * a2;
			s.q[i][0] = r0 * q0 - r1 * q1 - r2 * q2 - r3 * q3;
			s.q[i][1] = r0 * q1 + r1 * q0 + r2 * q3 - r3 * q2;
			s.q[i][2] = r0 * q2 - r1 * q3 + r2 * q0 + r3 * q1;
			s.q[i][3] = r0 * q3 + r1 * q2 - r2 * q1 + r3 * q0;
			assert(normalized(s.q[i]));
		}
	}

	// Calculate intra-ligand free energy.
	const size_t num_interacting_pairs = interacting_pairs.size();
	for (i = 0; i < num_interacting_pairs; ++i)
	{
		const interacting_pair& p = interacting_pairs[i];
		v0 = s.c[p.i2][0] - s.c[p.i1][0];
		v1 = s.c[p.i2][1] - s.c[p.i1][1];
		v2 = s.c[p.i2][2] - s.c[p.i1][2];
		vs = v0*v0 + v1*v1 + v2*v2;
		if (vs < scoring_function::cutoff_sqr)
		{
			o0 = p.p_offset + static_cast<size_t>(sf.ns * vs);
			e += sf.e[o0];
			dor = sf.d[o0];
			d0 = dor * v0;
			d1 = dor * v1;
			d2 = dor * v2;
			s.d[p.i1][0] -= d0;
			s.d[p.i1][1] -= d1;
			s.d[p.i1][2] -= d2;
			s.d[p.i2][0] += d0;
			s.d[p.i2][1] += d1;
			s.d[p.i2][2] += d2;
		}
	}

	// If the free energy is no better than the upper bound, refuse this conformation.
	if (e >= e_upper_bound) return false;

	// Store e from register into memory.
	s.e = e;

	// Calculate and aggregate the force and torque of BRANCH frames to their parent frame.
	for (i = 0; i < 3 * frames.size(); ++i)
	{
		s.f[i] = 0.0f;
		s.t[i] = 0.0f;
	}
	assert(k == frames.size());
	t = 6 + nt;
	while (true)
	{
		const frame& f = frames[--k];
		o0 = 3 * k;
		o1 = o0 + 1;
		o2 = o1 + 1;

		// Load variables from memory into registers.
		y0 = s.c[f.rotorYidx][0];
		y1 = s.c[f.rotorYidx][1];
		y2 = s.c[f.rotorYidx][2];
		f0 = s.f[o0];
		f1 = s.f[o1];
		f2 = s.f[o2];
		t0 = s.t[o0];
		t1 = s.t[o1];
		t2 = s.t[o2];
		for (i = f.beg; i < f.end; ++i)
		{
			d0 = s.d[i][0];
			d1 = s.d[i][1];
			d2 = s.d[i][2];

			// The derivatives with respect to the position, orientation, and torsions
			// would be the negative total force acting on the ligand,
			// the negative total torque, and the negative torque projections, respectively,
			// where the projections refer to the torque applied to the branch moved by the torsion,
			// projected on its rotation axis.
			f0 += d0;
			f1 += d1;
			f2 += d2;
			v0 = s.c[i][0] - y0;
			v1 = s.c[i][1] - y1;
			v2 = s.c[i][2] - y2;
			t0 += v1 * d2 - v2 * d1;
			t1 += v2 * d0 - v0 * d2;
			t2 += v0 * d1 - v1 * d0;
		}

		// Save the aggregated force and torque of ROOT frame to g.
		if (k == 0)
		{
			s.g[0] = f0;
			s.g[1] = f1;
			s.g[2] = f2;
			s.g[3] = t0;
			s.g[4] = t1;
			s.g[5] = t2;
			return true;
		}

		// Store variables from registers into memory.
		s.f[o0] = f0;
		s.f[o1] = f1;
		s.f[o2] = f2;
		s.t[o0] = t0;
		s.t[o1] = t1;
		s.t[o2] = t2;

		// Save the aggregated torque of active BRANCH frames to g.
		if (f.active)
		{
			s.g[--t] = t0 * s.a[o0] + t1 * s.a[o1] + t2 * s.a[o2]; // dot product
		}

		// Aggregate the force and torque of current frame to its parent frame.
		o0 = 3 * f.parent;
		o1 = o0 + 1;
		o2 = o1 + 1;
		s.f[o0] += f0;
		s.f[o1] += f1;
		s.f[o2] += f2;
		i = frames[f.parent].rotorYidx;
		v0 = y0 - s.c[i][0];
		v1 = y1 - s.c[i][1];
		v2 = y2 - s.c[i][2];
		s.t[o0] += t0 + v1 * f2 - v2 * f1;
		s.t[o1] += t1 + v2 * f0 - v0 * f2;
		s.t[o2] += t2 + v0 * f1 - v1 * f0;
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
	float q0, q1, q2, q3, qni, sum, alpha, pg1, pg2, po0, po1, po2, pon, hn, u, pq0, pq1, pq2, pq3, x1q0, x1q1, x1q2, x1q3, x2q0, x2q1, x2q2, x2q3, yhy, yp, ryp, pco;
	size_t g, i, j, o;
	mt19937_64 rng(seed);
	uniform_real_distribution<float> uniform_11(-1.0f, 1.0f);

	// Randomize s0.x.
	s0.x[0] = rec.center[0] + uniform_11(rng) * rec.size[0];
	s0.x[1] = rec.center[1] + uniform_11(rng) * rec.size[1];
	s0.x[2] = rec.center[2] + uniform_11(rng) * rec.size[2];
	q0 = uniform_11(rng);
	q1 = uniform_11(rng);
	q2 = uniform_11(rng);
	q3 = uniform_11(rng);
	qni = 1.0f / sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
	s0.x[3] = q0 * qni;
	s0.x[4] = q1 * qni;
	s0.x[5] = q2 * qni;
	s0.x[6] = q3 * qni;
	for (i = 0; i < nt; ++i)
	{
		s0.x[7 + i] = uniform_11(rng);
	}
	evaluate(s0, sf, rec, e_upper_bound);

	// Repeat for a number of generations.
	for (g = 0; g < num_generations; ++g)
	{
		// Mutate s0.x into s1.x
		s1.x[0] = s0.x[0] + uniform_11(rng);
		s1.x[1] = s0.x[1] + uniform_11(rng);
		s1.x[2] = s0.x[2] + uniform_11(rng);
		for (i = 3; i < nv + 1;  ++i)
		{
			s1.x[i] = s0.x[i];
		}
		evaluate(s1, sf, rec, e_upper_bound);

		// Initialize the inverse Hessian matrix to identity matrix.
		// An easier option that works fine in practice is to use a scalar multiple of the identity matrix,
		// where the scaling factor is chosen to be in the range of the eigenvalues of the true Hessian.
		// See N&R for a recipe to find this initializer.
		o = 0;
		for (j = 0; j < nv; ++j)
		{
			for (i = 0; i < j; ++i)
			{
				h[o++] = 0.0f;
			}
			h[o++] = 1.0f;
		}

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
				sum = 0.0f;
				for (j = 0; j < nv; ++j)
				{
					sum += h[mp(i, j)] * s1.g[j];
				}
				p[i] = -sum;
			}

			// Calculate pg = p*g = -h*g^2 < 0
			pg1 = 0;
			for (i = 0; i < nv; ++i)
			{
				pg1 += p[i] * s1.g[i];
			}

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
				po0 = p[3];
				po1 = p[4];
				po2 = p[5];
				pon = sqrt(po0*po0 + po1*po1 + po2*po2);
				hn = 0.5f * alpha * pon;
				u = sin(hn) / pon;
				pq0 = cos(hn);
				pq1 = u*po0;
				pq2 = u*po1;
				pq3 = u*po2;
				assert(fabs(pq0*pq0 + pq1*pq1 + pq2*pq2 + pq3*pq3 - 1.0f) < 1e-3f);
				x1q0 = s1.x[3];
				x1q1 = s1.x[4];
				x1q2 = s1.x[5];
				x1q3 = s1.x[6];
				assert(fabs(x1q0*x1q0 + x1q1*x1q1 + x1q2*x1q2 + x1q3*x1q3 - 1.0f) < 1e-3f);
				x2q0 = pq0 * x1q0 - pq1 * x1q1 - pq2 * x1q2 - pq3 * x1q3;
				x2q1 = pq0 * x1q1 + pq1 * x1q0 + pq2 * x1q3 - pq3 * x1q2;
				x2q2 = pq0 * x1q2 - pq1 * x1q3 + pq2 * x1q0 + pq3 * x1q1;
				x2q3 = pq0 * x1q3 + pq1 * x1q2 - pq2 * x1q1 + pq3 * x1q0;
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
					{
						pg2 += p[i] * s2.g[i];
					}
					if (pg2 >= 0.9f * pg1)
						break; // An appropriate alpha is found.
				}

				alpha *= 0.1f;
			}

			// If an appropriate alpha cannot be found, exit the BFGS loop.
			if (j == num_alphas) break;

			// Update Hessian matrix h.
			for (i = 0; i < nv; ++i) // Calculate y = g2 - g1.
			{
				y[i] = s2.g[i] - s1.g[i];
			}
			for (i = 0; i < nv; ++i) // Calculate mhy = -h * y.
			{
				sum = 0.0f;
				for (j = 0; j < nv; ++j)
				{
					sum += h[mp(i, j)] * y[j];
				}
				mhy[i] = -sum;
			}
			yhy = 0;
			for (i = 0; i < nv; ++i) // Calculate yhy = -y * mhy = -y * (-hy).
			{
				yhy -= y[i] * mhy[i];
			}
			yp = 0;
			for (i = 0; i < nv; ++i) // Calculate yp = y * p.
			{
				yp += y[i] * p[i];
			}
			ryp = 1 / yp;
			pco = ryp * (ryp * yhy + alpha);
			for (i = 0; i < nv; ++i)
			for (j = i; j < nv; ++j) // includes i
			{
				h[mr(i, j)] += ryp * (mhy[i] * p[j] + mhy[j] * p[i]) + pco * p[i] * p[j];
			}

			// Move to the next iteration.
			for (i = 0; i < nv + 1; ++i)
			{
				s1.x[i] = s2.x[i];
			}
			s1.e = s2.e;
			for (i = 0; i < nv; ++i)
			{
				s1.g[i] = s2.g[i];
			}
		}

		// Accept c1 according to Metropolis criteria.
		if (s1.e < s0.e)
		{
			for (i = 0; i < nv + 1; ++i)
			{
				s0.x[i] = s1.x[i];
			}
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
