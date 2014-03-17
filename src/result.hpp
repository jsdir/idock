#pragma once
#ifndef IDOCK_RESULT_HPP
#define IDOCK_RESULT_HPP

#include <boost/ptr_container/ptr_vector.hpp>
using boost::ptr_vector;

//! Represents a ligand conformation.
class conformation
{
public:
	array<double, 3> position; //!< Ligand origin coordinate.
	array<double, 4> orientation; //!< Ligand orientation.
	vector<double> torsions; //!< Ligand torsions.

	//! Constructs an initial conformation.
	explicit conformation(const size_t num_active_torsions) : position{}, orientation{1, 0, 0, 0}, torsions(num_active_torsions, 0) {}
};

//! Represents a transition from one conformation to another.
class change : public vector<double>
{
public:
	//! Constructs a zero change.
	explicit change(const size_t num_active_torsions) : vector<double>(6 + num_active_torsions, 0) {}
};

//! Represents a result found by BFGS local optimization for later clustering.
class result
{
public:
	double e; //!< Free energy.
	double f; //!< Inter-molecular free energy.
	double e_nd; //!< Normalized free energy, only for output purpose.
	vector<array<double, 3>> heavy_atoms; //!< Heavy atom coordinates.
	vector<array<double, 3>> hydrogens; //!< Hydrogen atom coordinates.

	//! Constructs a result from free energy e, force f, heavy atom coordinates and hydrogen atom coordinates.
	explicit result(const double e, const double f, vector<array<double, 3>>&& heavy_atoms_, vector<array<double, 3>>&& hydrogens_) : e(e), f(f), heavy_atoms(static_cast<vector<array<double, 3>>&&>(heavy_atoms_)), hydrogens(static_cast<vector<array<double, 3>>&&>(hydrogens_)) {}

	//! For sorting ptr_vector<result>.
	bool operator<(const result& r) const
	{
		return e < r.e;
	}
};

//! Clusters a result into an existing result set with a minimum RMSD requirement.
void add_to_result_container(ptr_vector<result>& results, result&& r, const double required_square_error);

#endif
