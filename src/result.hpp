#pragma once
#ifndef IDOCK_RESULT_HPP
#define IDOCK_RESULT_HPP

#include <boost/ptr_container/ptr_vector.hpp>
using boost::ptr_vector;

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
