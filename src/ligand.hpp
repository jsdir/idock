#pragma once
#ifndef IDOCK_LIGAND_HPP
#define IDOCK_LIGAND_HPP

#include <mutex>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/filesystem/fstream.hpp>
#include "atom.hpp"
#include "scoring_function.hpp"
#include "receptor.hpp"
#include "random_forest.hpp"
#include "mc_kernel.hpp"
using namespace boost::filesystem;

/// Represents a summary of docking results of a ligand.
class summary
{
public:
	const string stem;
	const vector<float> affinities;
	explicit summary(const string& stem, const vector<float>& affinities) : stem(stem), affinities(affinities) {}
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
	size_t childYidx; ///< The exclusive ending index to the heavy atoms of the current frame.
	bool active; ///< Indicates if the current frame is active.
	array<float, 3> yy; ///< Vector pointing from the origin of parent frame to the origin of current frame.
	array<float, 3> xy; ///< Normalized vector pointing from rotor X of parent frame to rotor Y of current frame.
	vector<size_t> branches; ///< Indexes to child branches.

	/// Constructs an active frame, and relates it to its parent frame.
	explicit frame(const size_t parent, const size_t rotorXsrn, const size_t rotorYsrn, const size_t rotorXidx, const size_t rotorYidx) : parent(parent), rotorXsrn(rotorXsrn), rotorYsrn(rotorYsrn), rotorXidx(rotorXidx), rotorYidx(rotorYidx), active(true) {}

	/// Outputs a BRANCH line in PDBQT format.
	void output(boost::filesystem::ofstream& ofs) const;
};

/// Represents a ligand.
class ligand
{
public:
	path p; ///< Path to the input ligand.
	vector<frame> frames; ///< ROOT and BRANCH frames.
	vector<atom> atoms; ///< Heavy atoms. Coordinates are relative to frame origin, which is the first atom by default.
	size_t nv; ///< Number of variables to optimize.
	size_t nf; ///< Number of frames.
	size_t na; ///< Number of atoms.
	size_t np; ///< Number of interacting pairs.
	vector<int> lig; ///< Encoded ligand content.

	/// Constructs a ligand by parsing a ligand file in pdbqt format.
	/// @exception parsing_error Thrown when an atom type is not recognized or an empty branch is detected.
	explicit ligand(const path p);

	int mc(const int tid, size_t& num_ligands, boost::ptr_vector<summary>& summaries, boost::ptr_vector<mc_kernel>& mc_kernels, const path& output_ligand_path, const size_t max_conformations, const size_t num_mc_tasks, const receptor& rec, const forest& f, mutex& m) const;

private:
	/// Represents a pair of interacting atoms that are separated by 3 consecutive covalent bonds.
	class interacting_pair
	{
	public:
		size_t i0; ///< Index of atom 0.
		size_t i1; ///< Index of atom 1.
		size_t p_offset; ///< Type pair index to the scoring function.
		interacting_pair(const size_t i0, const size_t i1, const size_t p_offset) : i0(i0), i1(i1), p_offset(p_offset) {}
	};

	vector<interacting_pair> interacting_pairs; ///< Non 1-4 interacting pairs.
};

#endif
