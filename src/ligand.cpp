#include <iomanip>
#include <random>
#include "utility.hpp"
#include "ligand.hpp"

void frame::output(boost::filesystem::ofstream& ofs) const
{
	ofs << "BRANCH"    << setw(4) << rotorXsrn << setw(4) << rotorYsrn << '\n';
}

ligand::ligand(const path& p) : num_active_torsions(0)
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
				++num_active_torsions;
			}

			// Set up bonds between rotorX and rotorY.
			bonds[f->rotorYidx].push_back(f->rotorXidx);
			bonds[f->rotorXidx].push_back(f->rotorYidx);

			// Dehydrophobicize rotorX and rotorY if necessary.
			atom& rotorY = atoms[f->rotorYidx];
			atom& rotorX = atoms[f->rotorXidx];
			if (rotorY.is_hetero() && !rotorX.is_hetero()) rotorX.dehydrophobicize();
			if (rotorX.is_hetero() && !rotorY.is_hetero()) rotorY.dehydrophobicize();

			// Calculate parent_rotorY_to_current_rotorY and parent_rotorX_to_current_rotorY.
			const frame& p = frames[f->parent];
			f->parent_rotorY_to_current_rotorY = rotorY.coord - atoms[p.rotorYidx].coord;
			f->parent_rotorX_to_current_rotorY = normalize(rotorY.coord - rotorX.coord);

			// Now the parent of the following frame is the parent of current frame.
			current = f->parent;

			// Update the pointer to the current frame.
			f = &frames[current];
		}
	}
	assert(current == 0); // current should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.
	assert(f == &frames.front()); // The frame pointer should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.
	assert(frames.size() >= 1);
	assert(frames.size() - 1 >= num_active_torsions);
	frames.back().end = atoms.size();
	num_variables = 6 + num_active_torsions;

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
			const vector<size_t>& i0_bonds = bonds[i];
			const size_t num_i0_bonds = i0_bonds.size();
			for (size_t i0 = 0; i0 < num_i0_bonds; ++i0)
			{
				const size_t b1 = i0_bonds[i0];
				if (find(neighbors.begin(), neighbors.end(), b1) == neighbors.end())
				{
					neighbors.push_back(b1);
				}
				const vector<size_t>& i1_bonds = bonds[b1];
				const size_t num_i1_bonds = i1_bonds.size();
				for (size_t i1 = 0; i1 < num_i1_bonds; ++i1)
				{
					const size_t b2 = i1_bonds[i1];
					if (find(neighbors.begin(), neighbors.end(), b2) == neighbors.end())
					{
						neighbors.push_back(b2);
					}
					const vector<size_t>& i2_bonds = bonds[b2];
					const size_t num_i2_bonds = i2_bonds.size();
					for (size_t i2 = 0; i2 < num_i2_bonds; ++i2)
					{
						const size_t b3 = i2_bonds[i2];
						if (find(neighbors.begin(), neighbors.end(), b3) == neighbors.end())
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

bool ligand::evaluate(const vector<float>& x, const scoring_function& sf, const receptor& rec, const float e_upper_bound, vector<array<float, 4>>& q, vector<array<float, 3>>& c, float& e, vector<float>& g) const
{
	vector<array<float, 3>> a(frames.size()); ///< Vector pointing from rotor Y to rotor X.
	vector<array<float, 3>> gf(frames.size(), zero3); ///< Aggregated derivatives of heavy atoms.
	vector<array<float, 3>> gt(frames.size(), zero3); /// Torque of the force.
	vector<array<float, 3>> d(atoms.size()); ///< Heavy atom derivatives.

	// Apply position and orientation to ROOT frame.
	const frame& root = frames.front();
	c.front()[0] = x[0];
	c.front()[1] = x[1];
	c.front()[2] = x[2];
	q.front()[0] = x[3];
	q.front()[1] = x[4];
	q.front()[2] = x[5];
	q.front()[3] = x[6];

	// Apply torsions to frames.
	for (size_t k = 0, t = 0; k < frames.size(); ++k)
	{
		const frame& f = frames[k];
		const array<float, 9> m = qtn4_to_mat3(q[k]);
		for (size_t i = f.beg + 1; i < f.end; ++i)
		{
			c[i] = c[f.rotorYidx] + m * atoms[i].coord;
		}
		for (const size_t i : f.branches)
		{
			const frame& b = frames[i];
			c[b.rotorYidx] = c[f.rotorYidx] + m * b.parent_rotorY_to_current_rotorY;

			// If the current BRANCH frame does not have an active torsion, skip it.
			if (!b.active)
			{
				assert(b.beg + 1 == b.end);
				assert(b.beg == b.rotorYidx);
				continue;
			}
			assert(normalized(b.parent_rotorX_to_current_rotorY));
			a[i] = m * b.parent_rotorX_to_current_rotorY;
			assert(normalized(a[i]));
			q[i] = vec4_to_qtn4(a[i], x[7 + t++]) * q[k];
			assert(normalized(q[i]));
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

	e = 0;
	for (size_t i = 0; i < atoms.size(); ++i)
	{
		if (!rec.within(c[i]))
		{
			e += 10;
			d[i][0] = 0;
			d[i][1] = 0;
			d[i][2] = 0;
			continue;
		}

		// Retrieve the grid map in need.
		const vector<float>& map = rec.maps[atoms[i].xs];
		assert(map.size());

		// Find the index of the current coordinates.
		const array<size_t, 3> index = rec.coordinate_to_index(c[i]);

		// Calculate the offsets to grid map and lookup the values.
		const size_t o000 = rec.num_probes[0] * (rec.num_probes[1] * index[2] + index[1]) + index[0];
		const size_t o100 = o000 + 1;
		const size_t o010 = o000 + rec.num_probes[0];
		const size_t o001 = o000 + rec.num_probes[0] * rec.num_probes[1];
		const float e000 = map[o000];
		const float e100 = map[o100];
		const float e010 = map[o010];
		const float e001 = map[o001];
		d[i][0] = (e100 - e000) * rec.granularity_inverse;
		d[i][1] = (e010 - e000) * rec.granularity_inverse;
		d[i][2] = (e001 - e000) * rec.granularity_inverse;

		e += e000; // Aggregate the energy.
	}

	// Calculate intra-ligand free energy.
	const size_t num_interacting_pairs = interacting_pairs.size();
	for (size_t i = 0; i < num_interacting_pairs; ++i)
	{
		const interacting_pair& p = interacting_pairs[i];
		const array<float, 3> r = c[p.i2] - c[p.i1];
		const float r2 = norm_sqr(r);
		if (r2 < scoring_function::cutoff_sqr)
		{
			const size_t o = p.p_offset + static_cast<size_t>(sf.ns * r2);
			e += sf.e[o];
			const array<float, 3> derivative = sf.d[o] * r;
			d[p.i1] -= derivative;
			d[p.i2] += derivative;
		}
	}

	// If the free energy is no better than the upper bound, refuse this conformation.
	if (e >= e_upper_bound) return false;

	// Calculate and aggregate the force and torque of BRANCH frames to their parent frame.
	for (size_t k = frames.size(), t = num_active_torsions; --k;)
	{
		const frame& f = frames[k];

		for (size_t i = f.beg; i < f.end; ++i)
		{
			// The derivatives with respect to the position, orientation, and torsions
			// would be the negative total force acting on the ligand,
			// the negative total torque, and the negative torque projections, respectively,
			// where the projections refer to the torque applied to the branch moved by the torsion,
			// projected on its rotation axis.
			gf[k] += d[i];
			gt[k] += (c[i] - c[f.rotorYidx]) * d[i];
		}

		// Aggregate the force and torque of current frame to its parent frame.
		gf[f.parent] += gf[k];
		gt[f.parent] += gt[k] + (c[f.rotorYidx] - c[frames[f.parent].rotorYidx]) * gf[k];

		// If the current BRANCH frame does not have an active torsion, skip it.
		if (!f.active) continue;

		// Save the torsion.
		g[6 + (--t)] = gt[k][0] * a[k][0] + gt[k][1] * a[k][1] + gt[k][2] * a[k][2]; // dot product
	}

	// Calculate and aggregate the force and torque of ROOT frame.
	for (size_t i = root.beg; i < root.end; ++i)
	{
		gf.front() += d[i];
		gt.front() += (c[i] - c.front()) * d[i];
	}

	// Save the aggregated force and torque to g.
	g[0] = gf.front()[0];
	g[1] = gf.front()[1];
	g[2] = gf.front()[2];
	g[3] = gt.front()[0];
	g[4] = gt.front()[1];
	g[5] = gt.front()[2];

	return true;
}

int ligand::bfgs(result& r, const scoring_function& sf, const receptor& rec, const size_t seed, const size_t num_generations) const
{
	// Define constants.
	const size_t num_alphas = 5; // Number of alpha values for determining step size in BFGS
	const float e_upper_bound = 40.0f * atoms.size(); // A conformation will be droped if its free energy is not better than e_upper_bound.

	// Declare variable.
	vector<float> x0(7 + num_active_torsions), x1(7 + num_active_torsions), x2(7 + num_active_torsions);
	vector<array<float, 4>> q0(frames.size()), q1(frames.size()), q2(frames.size());
	vector<array<float, 3>> c0(atoms.size()), c1(atoms.size()), c2(atoms.size());
	vector<float> g0(6 + num_active_torsions), g1(6 + num_active_torsions), g2(6 + num_active_torsions);
	vector<float> p(6 + num_active_torsions), y(6 + num_active_torsions), mhy(6 + num_active_torsions);
	vector<float> h(num_variables*(num_variables+1)>>1); // Symmetric triangular Hessian matrix.
	float e0, e1, e2, alpha, pg1, pg2, yhy, yp, ryp, pco;
	size_t g, i, j;
	mt19937_64 rng(seed);
	uniform_real_distribution<float> uniform_11(-1.0f, 1.0f);

	// Randomize conformation x0.
	x0[0] = rec.center[0] + uniform_11(rng) * rec.size[0];
	x0[1] = rec.center[1] + uniform_11(rng) * rec.size[1];
	x0[2] = rec.center[2] + uniform_11(rng) * rec.size[2];
	const array<float, 4> x0orientation = normalize(make_array(uniform_11(rng), uniform_11(rng), uniform_11(rng), uniform_11(rng)));
	assert(normalized(x0orientation));
	x0[3] = x0orientation[0];
	x0[4] = x0orientation[1];
	x0[5] = x0orientation[2];
	x0[6] = x0orientation[3];
	for (i = 0; i < num_active_torsions; ++i)
	{
		x0[7 + i] = uniform_11(rng);
	}
	evaluate(x0, sf, rec, e_upper_bound, q0, c0, e0, g0);
	r = result(q0, c0, e0);

	for (g = 0; g < num_generations; ++g)
	{
		// Make a copy, so the previous conformation is retained.
		x1 = x0;
		x1[0] += uniform_11(rng);
		x1[1] += uniform_11(rng);
		x1[2] += uniform_11(rng);
		evaluate(x1, sf, rec, e_upper_bound, q1, c1, e1, g1);

		// Initialize the inverse Hessian matrix to identity matrix.
		// An easier option that works fine in practice is to use a scalar multiple of the identity matrix,
		// where the scaling factor is chosen to be in the range of the eigenvalues of the true Hessian.
		// See N&R for a recipe to find this initializer.
		fill(h.begin(), h.end(), 0.0f);
		for (i = 0; i < num_variables; ++i)
			h[mr(i, i)] = 1.0f;

		// Given the mutated conformation c1, use BFGS to find a local minimum.
		// The conformation of the local minimum is saved to c2, and its derivative is saved to g2.
		// http://en.wikipedia.org/wiki/BFGS_method
		// http://en.wikipedia.org/wiki/Quasi-Newton_method
		// The loop breaks when an appropriate alpha cannot be found.
		while (true)
		{
			// Calculate p = -h*g, where p is for descent direction, h for Hessian, and g for gradient.
			for (i = 0; i < num_variables; ++i)
			{
				float sum = 0.0f;
				for (j = 0; j < num_variables; ++j)
					sum += h[mp(i, j)] * g1[j];
				p[i] = -sum;
			}

			// Calculate pg = p*g = -h*g^2 < 0
			pg1 = 0;
			for (i = 0; i < num_variables; ++i)
				pg1 += p[i] * g1[i];

			// Perform a line search to find an appropriate alpha.
			// Try different alpha values for num_alphas times.
			// alpha starts with 1, and shrinks to alpha_factor of itself iteration by iteration.
			alpha = 1.0;
			for (j = 0; j < num_alphas; ++j)
			{
				// Calculate c2 = c1 + ap.
				x2[0] = x1[0] + alpha * p[0];
				x2[1] = x1[1] + alpha * p[1];
				x2[2] = x1[2] + alpha * p[2];
				const array<float, 4> x1orientation = { x1[3], x1[4], x1[5], x1[6] };
				assert(normalized(x1orientation));
				const array<float, 4> x2orientation = vec3_to_qtn4(alpha * make_array(p[3], p[4], p[5])) * x1orientation;
				assert(normalized(x2orientation));
				x2[3] = x2orientation[0];
				x2[4] = x2orientation[1];
				x2[5] = x2orientation[2];
				x2[6] = x2orientation[3];
				for (i = 0; i < num_active_torsions; ++i)
				{
					x2[7 + i] = x1[7 + i] + alpha * p[6 + i];
				}

				// Evaluate c2, subject to Wolfe conditions http://en.wikipedia.org/wiki/Wolfe_conditions
				// 1) Armijo rule ensures that the step length alpha decreases f sufficiently.
				// 2) The curvature condition ensures that the slope has been reduced sufficiently.
				if (evaluate(x2, sf, rec, e1 + 0.0001f * alpha * pg1, q2, c2, e2, g2))
				{
					pg2 = 0;
					for (i = 0; i < num_variables; ++i)
						pg2 += p[i] * g2[i];
					if (pg2 >= 0.9f * pg1)
						break; // An appropriate alpha is found.
				}

				alpha *= 0.1f;
			}

			// If an appropriate alpha cannot be found, exit the BFGS loop.
			if (j == num_alphas) break;

			// Update Hessian matrix h.
			for (i = 0; i < num_variables; ++i) // Calculate y = g2 - g1.
				y[i] = g2[i] - g1[i];
			for (i = 0; i < num_variables; ++i) // Calculate mhy = -h * y.
			{
				float sum = 0.0f;
				for (j = 0; j < num_variables; ++j)
					sum += h[mp(i, j)] * y[j];
				mhy[i] = -sum;
			}
			yhy = 0;
			for (i = 0; i < num_variables; ++i) // Calculate yhy = -y * mhy = -y * (-hy).
				yhy -= y[i] * mhy[i];
			yp = 0;
			for (i = 0; i < num_variables; ++i) // Calculate yp = y * p.
				yp += y[i] * p[i];
			ryp = 1 / yp;
			pco = ryp * (ryp * yhy + alpha);
			for (i = 0; i < num_variables; ++i)
			for (j = i; j < num_variables; ++j) // includes i
			{
				h[mr(i, j)] += ryp * (mhy[i] * p[j] + mhy[j] * p[i]) + pco * p[i] * p[j];
			}

			// Move to the next iteration.
			x1 = x2;
			e1 = e2;
			g1 = g2;
		}

		// Accept c1 according to Metropolis criteria.
		if (e1 < e0)
		{
			r = result(q1, c1, e1);
			x0 = x1;
			e0 = e1;
		}
	}
	return 0;
}

void ligand::save(const path& output_ligand_path, const ptr_vector<result>& results, const vector<size_t>& representatives) const
{
	assert(representatives.size());
	assert(representatives.size() <= results.size());
	boost::filesystem::ofstream ofs(output_ligand_path);
	ofs.setf(ios::fixed, ios::floatfield);
	ofs << setprecision(3);
	for (size_t k = 0; k < representatives.size(); ++k)
	{
		const result& r = results[representatives[k]];
		ofs << "MODEL     " << setw(4) << (k + 1) << '\n'
			<< "REMARK     pKd:" << setw(7) << r.e << '\n';

		// Dump the ROOT frame.
		ofs << "ROOT\n";
		{
			const frame& f = frames.front();
			const array<float, 9> m = qtn4_to_mat3(r.q.front());
			for (size_t i = f.beg; i < f.end; ++i)
			{
				const atom& a = atoms[i];
				a.output(ofs, r.c[i]);
				if (a.hydrogens.empty()) continue;
				for (const atom& h : a.hydrogens)
				{
					h.output(ofs, r.c[f.rotorYidx] + m * h.coord);
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
				const array<float, 9> m = qtn4_to_mat3(f.active ? r.q[fn] : r.q[f.parent]);
				for (size_t i = f.beg; i < f.end; ++i)
				{
					const atom& a = atoms[i];
					a.output(ofs, r.c[i]);
					if (a.hydrogens.empty()) continue;
					for (const atom& h : a.hydrogens)
					{
						h.output(ofs, r.c[f.rotorYidx] + m * h.coord);
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
