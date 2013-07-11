#pragma once
#ifndef IDOCK_LIGAND_HPP
#define IDOCK_LIGAND_HPP

#include <boost/filesystem/fstream.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "atom.hpp"
#include "scoring_function.hpp"
#include "receptor.hpp"
using namespace boost::filesystem;
using boost::ptr_vector;

/// Represents a result found by BFGS local optimization for later clustering.
class solution
{
public:
	vector<float> x; ///< Conformation vector.
	vector<array<float, 3>> a; ///< Vector pointing from rotor Y to rotor X.
	vector<array<float, 4>> q; ///< Frame quaternions.
	vector<array<float, 3>> c; ///< Heavy atom coordinates.
	vector<array<float, 3>> d; ///< Heavy atom derivatives.
	vector<float> f; ///< Aggregated derivatives of heavy atoms.
	vector<float> t; /// Torque of the force.
	vector<float> g; ///< Gradient vector.
	float e; ///< Free energy.

	solution()
	{
	}

	void resize(const size_t nv, const size_t nf, const size_t na);

	/// For sorting ptr_vector<solution>.
	bool operator<(const solution& r) const
	{
		return e < r.e;
	}
};

/// Represents a ROOT or a BRANCH in PDBQT structure.
class frame
{
public:
	size_t parent; ///< Frame array index pointing to the parent of current frame. For ROOT frame, this field is not used.
	size_t rotorXsrn; ///< Serial atom number of the parent frame atom which forms a rotatable bond with the rotorY atom of current frame.
	size_t rotorYsrn; ///< Serial atom number of the current frame atom which forms a rotatable bond with the rotorX atom of parent frame.
	size_t rotorXidx; ///< Index pointing to the parent frame atom which forms a rotatable bond with the rotorY atom of current frame.
	size_t rotorYidx; ///< Index pointing to the current frame atom which forms a rotatable bond with the rotorX atom of parent frame.
	size_t beg; ///< The inclusive beginning index to the heavy atoms of the current frame.
	size_t end; ///< The exclusive ending index to the heavy atoms of the current frame.
	bool active; ///< Indicates if the current frame is active.
	array<float, 3> yy; ///< Vector pointing from the origin of parent frame to the origin of current frame.
	array<float, 3> xy; ///< Normalized vector pointing from rotor X of parent frame to rotor Y of current frame.
	vector<size_t> branches; ///< Indexes to child branches.

	/// Constructs an active frame, and relates it to its parent frame.
	explicit frame(const size_t parent, const size_t rotorXsrn, const size_t rotorYsrn, const size_t rotorXidx, const size_t beg) : parent(parent), rotorXsrn(rotorXsrn), rotorYsrn(rotorYsrn), rotorXidx(rotorXidx), beg(beg), active(true) {}

	/// Outputs a BRANCH line in PDBQT format.
	void output(boost::filesystem::ofstream& ofs) const;
};

/// Represents a ligand.
class ligand
{
public:
	vector<frame> frames; ///< ROOT and BRANCH frames.
	vector<atom> atoms; ///< Heavy atoms. Coordinates are relative to frame origin, which is the first atom by default.
	size_t nt; ///< Number of active torsions.
	size_t nv; ///< Number of variables to optimize.

	/// Constructs a ligand by parsing a ligand file in pdbqt format.
	/// @exception parsing_error Thrown when an atom type is not recognized or an empty branch is detected.
	ligand(const path& p);

	/// Evaluates free energy e, force f, and change g. Returns true if the conformation is accepted.
	bool evaluate(solution& s, const scoring_function& sf, const receptor& rec, const float e_upper_bound) const;

	/// Task for running Monte Carlo Simulated Annealing algorithm to find local minimums of the scoring function.
	int bfgs(solution& s, const scoring_function& sf, const receptor& rec, const size_t seed, const size_t num_generations) const;

	/// Writes a given number of conformations from a result container into a output ligand file in PDBQT format.
	void save(const path& output_ligand_path, const ptr_vector<solution>& solutions, const vector<size_t>& representatives) const;

private:
	/// Represents a pair of interacting atoms that are separated by 3 consecutive covalent bonds.
	class interacting_pair
	{
	public:
		size_t i1; ///< Index of atom 1.
		size_t i2; ///< Index of atom 2.
		size_t p_offset; ///< Type pair index to the scoring function.
		interacting_pair(const size_t i1, const size_t i2, const size_t p_offset) : i1(i1), i2(i2), p_offset(p_offset) {}
	};

	vector<interacting_pair> interacting_pairs; ///< Non 1-4 interacting pairs.
};

#endif
