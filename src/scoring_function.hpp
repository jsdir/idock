#pragma once
#ifndef IDOCK_SCORING_FUNCTION_HPP
#define IDOCK_SCORING_FUNCTION_HPP

#include <array>
#include <vector>
using namespace std;

//! Represents the scoring function used in idock.
class scoring_function
{
public:
	static const size_t n = 15; //!< Number of XScore atom types.
	static const size_t np = n*(n+1)>>1; //!< Number of XScore atom type pairs.
	static const size_t ns = 1024; //!< Number of samples in a unit distance.
	static const size_t cutoff = 8; //!< Atom type pair distance cutoff.
	static const size_t nr = ns*cutoff*cutoff+1; //!< Number of samples within the entire cutoff.
	static const size_t ne = nr*np; //!< Number of values to precalculate.
	static const float cutoff_sqr; //!< Cutoff square.

	//! Constructs an empty scoring function.
	explicit scoring_function();

	//! Aggregates the five term values evaluated at (t0, t1, r2).
	void score(float* const x, const size_t t0, const size_t t1, const float r2) const;

	//! Precalculates the scoring function values of sample points for the type combination of t0 and t1.
	void precalculate(const size_t t0, const size_t t1);

	//! Clears precalculated values.
	void clear();

	vector<float> e; //!< Scoring function values.
	vector<float> d; //!< Scoring function derivatives divided by distance.
private:
	static const array<float, n> vdw; //!< Van der Waals distances for XScore atom types.
	vector<float> rs; //!< Distance samples.
};

#endif
