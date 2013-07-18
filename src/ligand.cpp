#include <iomanip>
#include <random>
#include "utility.hpp"
#include "ligand.hpp"

void frame::output(boost::filesystem::ofstream& ofs) const
{
	ofs << "BRANCH"    << setw(4) << rotorXsrn << setw(4) << rotorYsrn << '\n';
}

ligand::ligand(const path& p) : nv(6)
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
				++nv;
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
	frames.back().end = na = atoms.size();
	nf = frames.size();
	assert(nf >= 1 + nv - 6);

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
	interacting_pairs.reserve(na * na);
	vector<size_t> neighbors;
	neighbors.reserve(10); // An atom typically consists of <= 10 neighbors.
	for (size_t k1 = 0; k1 < nf; ++k1)
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
			const size_t t1 = atoms[i].xs;
			for (size_t k2 = k1 + 1; k2 < nf; ++k2)
			{
				const frame& f2 = frames[k2];
				const frame& f3 = frames[f2.parent];
				for (size_t j = f2.beg; j < f2.end; ++j)
				{
					if (k1 == f2.parent && (i == f2.rotorXidx || j == f2.rotorYidx)) continue;
					if (k1 > 0 && f1.parent == f2.parent && i == f1.rotorYidx && j == f2.rotorYidx) continue;
					if (f2.parent > 0 && k1 == f3.parent && i == f3.rotorXidx && j == f2.rotorYidx) continue;
					if (find(neighbors.cbegin(), neighbors.cend(), j) != neighbors.cend()) continue;
					interacting_pairs.push_back(interacting_pair(i, j, scoring_function::nr * mp(t1, atoms[j].xs)));
				}
			}

			// Clear the current neighbor set for the next atom.
			neighbors.clear();
		}
	}
	np = interacting_pairs.size();

	lig.resize(11 * nf + nf - 1 + 4 * na + 3 * np);
	int* c = lig.data();
	for (const frame& f : frames) *c++ = f.active;
	for (const frame& f : frames) *c++ = f.beg;
	for (const frame& f : frames) *c++ = f.end;
	for (const frame& f : frames) *c++ = f.branches.size();
	for (const frame& f : frames) *c++ = f.parent;
	for (const frame& f : frames) *(float*)c++ = f.yy[0];
	for (const frame& f : frames) *(float*)c++ = f.yy[1];
	for (const frame& f : frames) *(float*)c++ = f.yy[2];
	for (const frame& f : frames) *(float*)c++ = f.xy[0];
	for (const frame& f : frames) *(float*)c++ = f.xy[1];
	for (const frame& f : frames) *(float*)c++ = f.xy[2];
	assert(c == lig.data() + 11 * nf);
	for (const frame& f : frames)
	{
		for (const size_t b : f.branches)
		{
			*c++ = b;
		}
	}
	assert(c == lig.data() + 11 * nf + nf - 1);
	for (const atom& a : atoms) *(float*)c++ = a.coord[0];
	for (const atom& a : atoms) *(float*)c++ = a.coord[1];
	for (const atom& a : atoms) *(float*)c++ = a.coord[2];
	for (const atom& a : atoms) *c++ = a.xs;
	assert(c == lig.data() + 11 * nf + nf - 1 + 4 * na);
	for (const interacting_pair& p : interacting_pairs) *c++ = p.i0;
	for (const interacting_pair& p : interacting_pairs) *c++ = p.i1;
	for (const interacting_pair& p : interacting_pairs) *c++ = p.p_offset;
	assert(c == lig.data() + 11 * nf + nf - 1 + 4 * na + 3 * np);
	assert(c == &lig.back() + 1);
}

bool ligand::evaluate(float* e, float* g, float* a, float* q, float* c, float* d, float* f, float* t, const float* x, const scoring_function& sf, const receptor& rec, const float eub, const unsigned int threadIdx, const unsigned int blockDim) const
{
	const int bd3 = 3 * blockDim;
	const int bd4 = 4 * blockDim;

	const int* act = lig.data();
	const int* beg = act + nf;
	const int* end = beg + nf;
	const int* nbr = end + nf;
	const int* prn = nbr + nf;
	const float* yy0 = (float*)(prn + nf);
	const float* yy1 = yy0 + nf;
	const float* yy2 = yy1 + nf;
	const float* xy0 = yy2 + nf;
	const float* xy1 = xy0 + nf;
	const float* xy2 = xy1 + nf;
	const int* brs = (int*)(xy2 + nf);
	const float* cr0 = (float*)(brs + nf - 1);
	const float* cr1 = cr0 + na;
	const float* cr2 = cr1 + na;
	const int* xst = (int*)(cr2 + na);
	const int* ip0 = xst + na;
	const int* ip1 = ip0 + np;
	const int* ipp = ip1 + np;
	assert(ipp + np == &lig.back() + 1);

	float y, y0, y1, y2, v0, v1, v2, c0, c1, c2, e000, e100, e010, e001, a0, a1, a2, ang, sng, r0, r1, r2, r3, vs, dor, f0, f1, f2, t0, t1, t2, d0, d1, d2;
	float q0, q1, q2, q3, q00, q01, q02, q03, q11, q12, q13, q22, q23, q33, m0, m1, m2, m3, m4, m5, m6, m7, m8;
	int i, j, k, b, w, i0, i1, i2, k0, k1, k2, z;

	// Apply position, orientation and torsions.
	c[i = threadIdx] = x[k = threadIdx];
	c[i += blockDim] = x[k += blockDim];
	c[i += blockDim] = x[k += blockDim];
	q[i = threadIdx] = x[k += blockDim];
	q[i += blockDim] = x[k += blockDim];
	q[i += blockDim] = x[k += blockDim];
	q[i += blockDim] = x[k += blockDim];
	y = 0.0f;
	for (k = 0, b = 0, w = 6 * blockDim + threadIdx; k < nf; ++k)
	{
		// Load rotorY from memory into registers.
		y0 = c[i0 = beg[k] * bd3 + threadIdx];
		y1 = c[i0 += blockDim];
		y2 = c[i0 += blockDim];

		// Translate orientation of active frames from quaternion into 3x3 matrix.
		if (act[k])
		{
			q0 = q[k0 = k * bd4 + threadIdx];
			q1 = q[k0 += blockDim];
			q2 = q[k0 += blockDim];
			q3 = q[k0 += blockDim];
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

		// Evaluate c and d of frame atoms. Aggregate e into y.
		for (i = beg[k], z = end[k]; i < z; ++i)
		{
			i0 = i * bd3 + threadIdx;
			i1 = i0 + blockDim;
			i2 = i1 + blockDim;

			// The first atom of a frame is assumed to be its rotor Y.
			if (i == beg[k])
			{
				c0 = y0;
				c1 = y1;
				c2 = y2;
			}
			else
			{
				// Calculate coordinate from transformation matrix and offset.
				v0 = cr0[i];
				v1 = cr1[i];
				v2 = cr2[i];
				c0 = y0 + m0 * v0 + m1 * v1 + m2 * v2;
				c1 = y1 + m3 * v0 + m4 * v1 + m5 * v2;
				c2 = y2 + m6 * v0 + m7 * v1 + m8 * v2;

				// Store coordinate from registers into memory.
				c[i0] = c0;
				c[i1] = c1;
				c[i2] = c2;
			}

			// Penalize out-of-box case.
			if (c0 < rec.corner0[0] || rec.corner1[0] <= c0 || c1 < rec.corner0[1] || rec.corner1[1] <= c1 || c2 < rec.corner0[2] || rec.corner1[2] <= c2)
			{
				y += 10.0f;
				d[i0] = 0.0f;
				d[i1] = 0.0f;
				d[i2] = 0.0f;
				continue;
			}

			// Find the index of the current coordinate
			k0 = static_cast<size_t>((c0 - rec.corner0[0]) * rec.granularity_inverse);
			k1 = static_cast<size_t>((c1 - rec.corner0[1]) * rec.granularity_inverse);
			k2 = static_cast<size_t>((c2 - rec.corner0[2]) * rec.granularity_inverse);
			assert(k0 + 1 < rec.num_probes[0]);
			assert(k1 + 1 < rec.num_probes[1]);
			assert(k2 + 1 < rec.num_probes[2]);
			k0 = rec.num_probes[0] * (rec.num_probes[1] * k2 + k1) + k0;

			// Retrieve the grid map and lookup the value
			const vector<float>& map = rec.maps[xst[i]];
			assert(map.size());
			e000 = map[k0];
			e100 = map[k0 + 1];
			e010 = map[k0 + rec.num_probes[0]];
			e001 = map[k0 + rec.num_probes[0] * rec.num_probes[1]];
			y += e000;
			d[i0] = (e100 - e000) * rec.granularity_inverse;
			d[i1] = (e010 - e000) * rec.granularity_inverse;
			d[i2] = (e001 - e000) * rec.granularity_inverse;
		}
		for (j = 0, z = nbr[k]; j < z; ++j)
		{
			i = brs[b++];
			i0 = beg[i] * bd3 + threadIdx;
			i1 = i0 + blockDim;
			i2 = i1 + blockDim;
			c[i0] = y0 + m0 * yy0[i] + m1 * yy1[i] + m2 * yy2[i];
			c[i1] = y1 + m3 * yy0[i] + m4 * yy1[i] + m5 * yy2[i];
			c[i2] = y2 + m6 * yy0[i] + m7 * yy1[i] + m8 * yy2[i];

			// Skip inactive BRANCH frame
			if (!act[i]) continue;

			// Update a of BRANCH frame
			a0 = m0 * xy0[i] + m1 * xy1[i] + m2 * xy2[i];
			a1 = m3 * xy0[i] + m4 * xy1[i] + m5 * xy2[i];
			a2 = m6 * xy0[i] + m7 * xy1[i] + m8 * xy2[i];
			assert(fabs(a0*a0 + a1*a1 + a2*a2 - 1.0f) < 1e-3f);
			a[k0 = i * bd3 + threadIdx] = a0;
			a[k0 += blockDim] = a1;
			a[k0 += blockDim] = a2;

			// Update q of BRANCH frame
			ang = x[w += blockDim] * 0.5f;
			sng = sinf(ang);
			r0 = cosf(ang);
//			sincosf(ang, &sng, &r0);
			r1 = sng * a0;
			r2 = sng * a1;
			r3 = sng * a2;
			q00 = r0 * q0 - r1 * q1 - r2 * q2 - r3 * q3;
			q01 = r0 * q1 + r1 * q0 + r2 * q3 - r3 * q2;
			q02 = r0 * q2 - r1 * q3 + r2 * q0 + r3 * q1;
			q03 = r0 * q3 + r1 * q2 - r2 * q1 + r3 * q0;
			assert(fabs(q00*q00 + q01*q01 + q02*q02 + q03*q03 - 1.0f) < 1e-3f);
			q[k0 = i * bd4 + threadIdx] = q00;
			q[k0 += blockDim] = q01;
			q[k0 += blockDim] = q02;
			q[k0 += blockDim] = q03;
		}
	}
	assert(b == nf - 1);
	assert(w == nv * blockDim + threadIdx);
	assert(k == nf);

	// Calculate intra-ligand free energy.
	for (i = 0; i < np; ++i)
	{
		i0 = ip0[i] * bd3 + threadIdx;
		i1 = i0 + blockDim;
		i2 = i1 + blockDim;
		k0 = ip1[i] * bd3 + threadIdx;
		k1 = k0 + blockDim;
		k2 = k1 + blockDim;
		v0 = c[k0] - c[i0];
		v1 = c[k1] - c[i1];
		v2 = c[k2] - c[i2];
		vs = v0*v0 + v1*v1 + v2*v2;
		if (vs < scoring_function::cutoff_sqr)
		{
			j = ipp[i] + static_cast<size_t>(sf.ns * vs);
			y += sf.e[j];
			dor = sf.d[j];
			d0 = dor * v0;
			d1 = dor * v1;
			d2 = dor * v2;
			d[i0] -= d0;
			d[i1] -= d1;
			d[i2] -= d2;
			d[k0] += d0;
			d[k1] += d1;
			d[k2] += d2;
		}
	}

	// If the free energy is no better than the upper bound, refuse this conformation.
	if (y >= eub) return false;

	// Store e from register into memory.
	e[threadIdx] = y;

	// Calculate and aggregate the force and torque of BRANCH frames to their parent frame.
	f[k0 = threadIdx] = 0.0f;
	t[k0] = 0.0f;
	for (i = 1, z = 3 * nf; i < z; ++i)
	{
		f[k0 += blockDim] = 0.0f;
		t[k0] = 0.0f;
	}
	assert(w == nv * blockDim + threadIdx);
	assert(k == nf);
	while (true)
	{
		--k;

		// Load f, t and rotorY from memory into register
		k0 = k * bd3 + threadIdx;
		k1 = k0 + blockDim;
		k2 = k1 + blockDim;
		f0 = f[k0];
		f1 = f[k1];
		f2 = f[k2];
		t0 = t[k0];
		t1 = t[k1];
		t2 = t[k2];
		y0 = c[i0 = beg[k] * bd3 + threadIdx];
		y1 = c[i0 += blockDim];
		y2 = c[i0 += blockDim];

		// Aggregate frame atoms.
		for (i = beg[k], z = end[k]; i < z; ++i)
		{
			i0 = i * bd3 + threadIdx;
			i1 = i0 + blockDim;
			i2 = i1 + blockDim;
			d0 = d[i0];
			d1 = d[i1];
			d2 = d[i2];

			// The derivatives with respect to the position, orientation, and torsions
			// would be the negative total force acting on the ligand,
			// the negative total torque, and the negative torque projections, respectively,
			// where the projections refer to the torque applied to the branch moved by the torsion,
			// projected on its rotation axi
			f0 += d0;
			f1 += d1;
			f2 += d2;
			if (i == beg[k]) continue;

			v0 = c[i0] - y0;
			v1 = c[i1] - y1;
			v2 = c[i2] - y2;
			t0 += v1 * d2 - v2 * d1;
			t1 += v2 * d0 - v0 * d2;
			t2 += v0 * d1 - v1 * d0;
		}

		// Save the aggregated force and torque of ROOT frame to g.
		if (k == 0)
		{
			g[i0 = threadIdx] = f0;
			g[i0 += blockDim] = f1;
			g[i0 += blockDim] = f2;
			g[i0 += blockDim] = t0;
			g[i0 += blockDim] = t1;
			g[i0 += blockDim] = t2;
			return true;
		}

		// Save the aggregated torque of active BRANCH frames to g.
		if (act[k])
		{
			g[w -= blockDim] = t0 * a[k0] + t1 * a[k1] + t2 * a[k2]; // dot product
		}

		// Aggregate the force and torque of current frame to its parent frame.
		k0 = prn[k] * bd3 + threadIdx;
		k1 = k0 + blockDim;
		k2 = k1 + blockDim;
		f[k0] += f0;
		f[k1] += f1;
		f[k2] += f2;
		v0 = y0 - c[i0 = beg[prn[k]] * bd3 + threadIdx];
		v1 = y1 - c[i0 += blockDim];
		v2 = y2 - c[i0 += blockDim];
		t[k0] += t0 + v1 * f2 - v2 * f1;
		t[k1] += t1 + v2 * f0 - v0 * f2;
		t[k2] += t2 + v0 * f1 - v1 * f0;
	}
	assert(w == 6 * blockDim + threadIdx);
}

int ligand::bfgs(float* s0e, const scoring_function& sf, const receptor& rec, const size_t seed, const size_t ng, const unsigned int threadIdx, const unsigned int blockDim) const
{
	const int nls = 5; // Number of line search trials for determining step size in BFGS
	const float eub = 40.0f * na; // A conformation will be droped if its free energy is not better than e_upper_bound.
	float* s0x = s0e + blockDim;
	float* s0g = s0x + (nv + 1) * blockDim;
	float* s0a = s0g + nv * blockDim;
	float* s0q = s0a + 3 * nf * blockDim;
	float* s0c = s0q + 4 * nf * blockDim;
	float* s0d = s0c + 3 * na * blockDim;
	float* s0f = s0d + 3 * na * blockDim;
	float* s0t = s0f + 3 * nf * blockDim;
	float* s1e = s0t + 3 * nf * blockDim;
	float* s1x = s1e + blockDim;
	float* s1g = s1x + (nv + 1) * blockDim;
	float* s1a = s1g + nv * blockDim;
	float* s1q = s1a + 3 * nf * blockDim;
	float* s1c = s1q + 4 * nf * blockDim;
	float* s1d = s1c + 3 * na * blockDim;
	float* s1f = s1d + 3 * na * blockDim;
	float* s1t = s1f + 3 * nf * blockDim;
	float* s2e = s1t + 3 * nf * blockDim;
	float* s2x = s2e + blockDim;
	float* s2g = s2x + (nv + 1) * blockDim;
	float* s2a = s2g + nv * blockDim;
	float* s2q = s2a + 3 * nf * blockDim;
	float* s2c = s2q + 4 * nf * blockDim;
	float* s2d = s2c + 3 * na * blockDim;
	float* s2f = s2d + 3 * na * blockDim;
	float* s2t = s2f + 3 * nf * blockDim;
	float* bfh = s2t + 3 * nf * blockDim;
	float* bfp = bfh + (nv*(nv+1)>>1) * blockDim;
	float* bfy = bfp + nv * blockDim;
	float* bfm = bfy + nv * blockDim;
	float rd0, rd1, rd2, rd3, rst;
	float sum, pg1, pga, pgc, alp, pg2, pr0, pr1, pr2, nrm, ang, sng, pq0, pq1, pq2, pq3, s1xq0, s1xq1, s1xq2, s1xq3, s2xq0, s2xq1, s2xq2, s2xq3, bpi;
	float yhy, yps, ryp, pco, bpj, bmj, ppj;
	int g, i, j, o0, o1, o2;
	mt19937_64 rng(seed);
	uniform_real_distribution<float> uniform_11(-1.0f, 1.0f);

	// Randomize s0x.
	rd0 = uniform_11(rng);
	s0x[o0 = threadIdx] = 0.5f * ((1 + rd0) * rec.corner1[0] + (1 - rd0) * rec.corner0[0]);
	rd0 = uniform_11(rng);
	s0x[o0 += blockDim] = 0.5f * ((1 + rd0) * rec.corner1[1] + (1 - rd0) * rec.corner0[1]);
	rd0 = uniform_11(rng);
	s0x[o0 += blockDim] = 0.5f * ((1 + rd0) * rec.corner1[2] + (1 - rd0) * rec.corner0[2]);
	rd0 = uniform_11(rng);
	rd1 = uniform_11(rng);
	rd2 = uniform_11(rng);
	rd3 = uniform_11(rng);
	rst = 1.0f / sqrt(rd0*rd0 + rd1*rd1 + rd2*rd2 + rd3*rd3);
//	rst = rsqrtf(rd0*rd0 + rd1*rd1 + rd2*rd2 + rd3*rd3);
	s0x[o0 += blockDim] = rd0 * rst;
	s0x[o0 += blockDim] = rd1 * rst;
	s0x[o0 += blockDim] = rd2 * rst;
	s0x[o0 += blockDim] = rd3 * rst;
	for (i = 6; i < nv; ++i)
	{
		s0x[o0 += blockDim] = uniform_11(rng);
	}
/*
	s0x[o0 = threadIdx] =  49.799f;
	s0x[o0 += blockDim] = -31.025f;
	s0x[o0 += blockDim] =  35.312f;
	s0x[o0 += blockDim] = 1.0f;
	s0x[o0 += blockDim] = 0.0f;
	s0x[o0 += blockDim] = 0.0f;
	s0x[o0 += blockDim] = 0.0f;
	for (i = 6; i < nv; ++i)
	{
		s0x[o0 += blockDim] = 0.0f;
	}
*/
	evaluate(s0e, s0g, s0a, s0q, s0c, s0d, s0f, s0t, s0x, sf, rec, eub, threadIdx, blockDim);

	// Repeat for a number of generations.
	for (g = 0; g < ng; ++g)
	{
		// Mutate s0x into s1x
		o0 = threadIdx;
		s1x[o0] = s0x[o0] + uniform_11(rng);
		o0 += blockDim;
		s1x[o0] = s0x[o0] + uniform_11(rng);
		o0 += blockDim;
		s1x[o0] = s0x[o0] + uniform_11(rng);
//		for (i = 3; i < nv + 1; ++i)
		for (i = 2 - nv; i < 0; ++i)
		{
			o0 += blockDim;
			s1x[o0] = s0x[o0];
		}
		evaluate(s1e, s1g, s1a, s1q, s1c, s1d, s1f, s1t, s1x, sf, rec, eub, threadIdx, blockDim);

		// Initialize the inverse Hessian matrix to identity matrix.
		// An easier option that works fine in practice is to use a scalar multiple of the identity matrix,
		// where the scaling factor is chosen to be in the range of the eigenvalues of the true Hessian.
		// See N&R for a recipe to find this initializer.
		bfh[o0 = threadIdx] = 1.0f;
		for (j = 1; j < nv; ++j)
		{
			for (i = 0; i < j; ++i)
			{
				bfh[o0 += blockDim] = 0.0f;
			}
			bfh[o0 += blockDim] = 1.0f;
		}

		// Use BFGS to optimize the mutated conformation s1x into local optimum s2x.
		// http://en.wikipedia.org/wiki/BFGS_method
		// http://en.wikipedia.org/wiki/Quasi-Newton_method
		// The loop breaks when no appropriate alpha can be found.
		while (true)
		{
			// Calculate p = -h * g, where p is for descent direction, h for Hessian, and g for gradient.
			sum = bfh[o1 = threadIdx] * s1g[o0 = threadIdx];
			for (i = 1; i < nv; ++i)
			{
				sum += bfh[o1 += i * blockDim] * s1g[o0 += blockDim];
			}
			bfp[o2 = threadIdx] = -sum;
			for (j = 1; j < nv; ++j)
			{
				sum = bfh[o1 = (j*(j+1)>>1) * blockDim + threadIdx] * s1g[o0 = threadIdx];
				for (i = 1; i < nv; ++i)
				{
					sum += bfh[o1 += i > j ? i * blockDim : blockDim] * s1g[o0 += blockDim];
				}
				bfp[o2 += blockDim] = -sum;
			}

			// Calculate pg = p * g = -h * g^2 < 0
			o0 = threadIdx;
			pg1 = bfp[o0] * s1g[o0];
			for (i = 1; i < nv; ++i)
			{
				o0 += blockDim;
				pg1 += bfp[o0] * s1g[o0];
			}
			pga = 0.0001f * pg1;
			pgc = 0.9f * pg1;

			// Perform a line search to find an appropriate alpha.
			// Try different alpha values for nls times.
			// alpha starts with 1, and shrinks to 0.1 of itself iteration by iteration.
			alp = 1.0f;
			for (j = 0; j < nls; ++j)
			{
				// Calculate x2 = x1 + a * p.
				o0 = threadIdx;
				s2x[o0] = s1x[o0] + alp * bfp[o0];
				o0 += blockDim;
				s2x[o0] = s1x[o0] + alp * bfp[o0];
				o0 += blockDim;
				s2x[o0] = s1x[o0] + alp * bfp[o0];
				o0 += blockDim;
				s1xq0 = s1x[o0];
				pr0 = bfp[o0];
				o0 += blockDim;
				s1xq1 = s1x[o0];
				pr1 = bfp[o0];
				o0 += blockDim;
				s1xq2 = s1x[o0];
				pr2 = bfp[o0];
				o0 += blockDim;
				s1xq3 = s1x[o0];
				assert(fabs(s1xq0*s1xq0 + s1xq1*s1xq1 + s1xq2*s1xq2 + s1xq3*s1xq3 - 1.0f) < 1e-3f);
				nrm = sqrt(pr0*pr0 + pr1*pr1 + pr2*pr2);
				ang = 0.5f * alp * nrm;
				sng = sinf(ang) / nrm;
				pq0 = cosf(ang);
//				sincosf(ang, &sng, &pq0);
//				sng /= nrm;
				pq1 = sng * pr0;
				pq2 = sng * pr1;
				pq3 = sng * pr2;
				assert(fabs(pq0*pq0 + pq1*pq1 + pq2*pq2 + pq3*pq3 - 1.0f) < 1e-3f);
				s2xq0 = pq0 * s1xq0 - pq1 * s1xq1 - pq2 * s1xq2 - pq3 * s1xq3;
				s2xq1 = pq0 * s1xq1 + pq1 * s1xq0 + pq2 * s1xq3 - pq3 * s1xq2;
				s2xq2 = pq0 * s1xq2 - pq1 * s1xq3 + pq2 * s1xq0 + pq3 * s1xq1;
				s2xq3 = pq0 * s1xq3 + pq1 * s1xq2 - pq2 * s1xq1 + pq3 * s1xq0;
				assert(fabs(s2xq0*s2xq0 + s2xq1*s2xq1 + s2xq2*s2xq2 + s2xq3*s2xq3 - 1.0f) < 1e-3f);
				s2x[o0 -= 3 * blockDim] = s2xq0;
				s2x[o0 += blockDim] = s2xq1;
				s2x[o0 += blockDim] = s2xq2;
				s2x[o0 += blockDim] = s2xq3;
				for (i = 6; i < nv; ++i)
				{
					bpi = bfp[o0];
					o0 += blockDim;
					s2x[o0] = s1x[o0] + alp * bpi;
				}

				// Evaluate x2, subject to Wolfe conditions http://en.wikipedia.org/wiki/Wolfe_conditions
				// 1) Armijo rule ensures that the step length alpha decreases f sufficiently.
				// 2) The curvature condition ensures that the slope has been reduced sufficiently.
				if (evaluate(s2e, s2g, s2a, s2q, s2c, s2d, s2f, s2t, s2x, sf, rec, s1e[threadIdx] + alp * pga, threadIdx, blockDim))
				{
					o0 = threadIdx;
					pg2 = bfp[o0] * s2g[o0];
					for (i = 1; i < nv; ++i)
					{
						o0 += blockDim;
						pg2 += bfp[o0] * s2g[o0];
					}
					if (pg2 >= pgc) break;
				}

				alp *= 0.1f;
			}

			// If no appropriate alpha can be found, exit the BFGS loop.
			if (j == nls) break;

			// Calculate y = g2 - g1.
			o0 = threadIdx;
			bfy[o0] = s2g[o0] - s1g[o0];
			for (i = 1; i < nv; ++i)
			{
				o0 += blockDim;
				bfy[o0] = s2g[o0] - s1g[o0];
			}

			// Calculate m = -h * y.
			sum = bfh[o1 = threadIdx] * bfy[o0 = threadIdx];
			for (i = 1; i < nv; ++i)
			{
				sum += bfh[o1 += i * blockDim] * bfy[o0 += blockDim];
			}
			bfm[o2 = threadIdx] = -sum;
			for (j = 1; j < nv; ++j)
			{
				sum = bfh[o1 = (j*(j+1)>>1) * blockDim + threadIdx] * bfy[o0 = threadIdx];
				for (i = 1; i < nv; ++i)
				{
					sum += bfh[o1 += i > j ? i * blockDim : blockDim] * bfy[o0 += blockDim];
				}
				bfm[o2 += blockDim] = -sum;
			}

			// Calculate yhy = -y * m = -y * (-h * y) = y * h * y.
			o0 = threadIdx;
			yhy = -bfy[o0] * bfm[o0];
			for (i = 1; i < nv; ++i)
			{
				o0 += blockDim;
				yhy -= bfy[o0] * bfm[o0];
			}

			// Calculate yps = y * p.
			o0 = threadIdx;
			yps = bfy[o0] * bfp[o0];
			for (i = 1; i < nv; ++i)
			{
				o0 += blockDim;
				yps += bfy[o0] * bfp[o0];
			}

			// Update Hessian matrix h.
			ryp = 1 / yps;
			pco = ryp * (ryp * yhy + alp);
			o2 = threadIdx;
			for (j = 0; j < nv; ++j)
			{
				bpj = bfp[o2];
				bmj = bfm[o2];
				ppj = pco * bpj;
				bfh[o1 = (j*(j+3)>>1) * blockDim + threadIdx] += (ryp * 2 * bmj + ppj) * bpj;
				for (i = j + 1; i < nv; ++i)
				{
					o0 = i * blockDim + threadIdx;
					bpi = bfp[o0];
					bfh[o1 += i * blockDim] += ryp * (bmj * bpi + bfm[o0] * bpj) + ppj * bpi;
				}
				o2 += blockDim;
			}

			// Move to the next iteration, i.e. e1 = e2, x1 = x2, g1 = g2.
			o0 = threadIdx;
			s1e[o0] = s2e[o0];
//			for (i = 1; i < 2 * (nv + 1); ++i)
			for (i = -1 - 2 * nv; i < 0; ++i)
			{
				o0 += blockDim;
				s1e[o0] = s2e[o0];
			}
		}

		// Accept x1 according to Metropolis criteria.
		if (s1e[threadIdx] < s0e[threadIdx])
		{
			o0 = threadIdx;
			s0e[o0] = s1e[o0];
//			for (i = 1; i < nv + 2; ++i)
			for (i = -1 - nv; i < 0; ++i)
			{
				o0 += blockDim;
				s0e[o0] = s1e[o0];
			}
		}
	}
	return 0;
}

void ligand::recover(solution& s) const
{
	// Apply position and orientation to ROOT frame.
	s.q.resize(nf);
	s.c.resize(na);
	s.c[0][0] = s.x[0];
	s.c[0][1] = s.x[1];
	s.c[0][2] = s.x[2];
	s.q[0][0] = s.x[3];
	s.q[0][1] = s.x[4];
	s.q[0][2] = s.x[5];
	s.q[0][3] = s.x[6];

	// Apply torsions to frames.
	for (size_t k = 0, w = 7; k < nf; ++k)
	{
		const frame& f = frames[k];
		if (!f.active) continue;
		const array<float, 9> m = qtn4_to_mat3(s.q[k]);
		for (size_t i = f.beg + 1; i < f.end; ++i)
		{
			s.c[i] = s.c[f.beg] + m * atoms[i].coord;
		}
		for (const size_t i : f.branches)
		{
			const frame& b = frames[i];
			s.c[b.beg] = s.c[f.beg] + m * b.yy;

			// Skip inactive BRANCH frame
			if (!b.active) continue;

			const array<float, 3> a = m * b.xy;
			assert(normalized(a));
			s.q[i] = vec4_to_qtn4(a, s.x[w++]) * s.q[k];
			assert(normalized(s.q[i]));
		}
	}
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

		// Dump the ROOT frame.
		ofs << "ROOT\n";
		{
			const frame& f = frames.front();
			const array<float, 9> m = qtn4_to_mat3(s.q[0]);
			for (size_t i = f.beg; i < f.end; ++i)
			{
				const atom& a = atoms[i];
				a.output(ofs, s.c[i]);
				for (const atom& h : a.hydrogens)
				{
					h.output(ofs, s.c[f.beg] + m * h.coord);
				}
			}
		}
		ofs << "ENDROOT\n";

		// Dump the BRANCH frames.
		vector<bool> dumped(nf); // dump_branches[0] is dummy. The ROOT frame has been dumped.
		vector<size_t> stack; // Stack to track the depth-first traversal sequence of frames in order to avoid recursion.
		stack.reserve(nf - 1); // The ROOT frame is excluded.
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
				const array<float, 9> m = qtn4_to_mat3(s.q[f.active ? fn : f.parent]);
				for (size_t i = f.beg; i < f.end; ++i)
				{
					const atom& a = atoms[i];
					a.output(ofs, s.c[i]);
					for (const atom& h : a.hydrogens)
					{
						h.output(ofs, s.c[f.beg] + m * h.coord);
					}
				}
				dumped[fn] = true;
				for (auto i = f.branches.rbegin(); i < f.branches.rend(); ++i)
				{
					stack.push_back(*i);
				}
			}
		}
		ofs << "TORSDOF " << nf - 1 << '\n';
	}
}
