#pragma once
#ifndef IDOCK_SCORING_FUNCTION_HPP
#define IDOCK_SCORING_FUNCTION_HPP

#include "atom.hpp"
#include "matrix.hpp"

//! Represents a pair of scoring function value and dor at a specific combination of (t0, t1, r).
class scoring_function_element
{
public:
	double e; //!< Scoring function value.
	double dor; //!< Scoring function derivative over r.
};

//! Represents a scoring function.
class scoring_function : public triangular_matrix<vector<scoring_function_element>>
{
public:
	static const size_t n = 15; //!< Number of XScore atom types.
	static const double Cutoff; //!< Cutoff of a scoring function.
	static const double Cutoff_Sqr; //!< Square of Cutoff.
	static const double Factor; //!< Scaling factor for r, i.e. distance between two atoms.
	static const double Factor_Inverse; //!< 1 / Factor.
	static const size_t Num_Samples; //!< Number of sampling points within [0, Cutoff].

	//! Constructs an empty scoring function.
	explicit scoring_function();

	//! Returns the score between two atoms of XScore atom types t0 and t1 and distance r.
	static double score(const size_t t0, const size_t t1, const double r);

	static void score(double* const v, const size_t t0, const size_t t1, const double r2);

	//! Precalculates the scoring function values of sample points for the type combination of t0 and t1.
	void precalculate(const size_t t0, const size_t t1);

	//! Evaluates the scoring function given (t0, t1, r2).
	scoring_function_element evaluate(const size_t type_pair_index, const double r2) const;

	//! Clears precalculated values.
	void clear();

private:
	static const array<double, n> vdw; //!< Van der Waals distances for XScore atom types.
	vector<double> rs; //!< Distance samples.
};

#endif
