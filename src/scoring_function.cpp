#include <cmath>
#include "scoring_function.hpp"

const double scoring_function::Cutoff = static_cast<double>(8);
const double scoring_function::Cutoff_Sqr = Cutoff * Cutoff;
const double scoring_function::Factor = static_cast<double>(256);
const double scoring_function::Factor_Inverse = 1 / Factor;
const size_t scoring_function::Num_Samples = static_cast<size_t>(Factor * Cutoff_Sqr) + 1;

const array<float, scoring_function::n> scoring_function::vdw =
{
	1.9, //   C_H
	1.9, //   C_P
	1.8, //   N_P
	1.8, //   N_D
	1.8, //   N_A
	1.8, //   N_DA
	1.7, //   O_A
	1.7, //   O_DA
	2.0, //   S_P
	2.1, //   P_P
	1.5, //   F_H
	1.8, //  Cl_H
	2.0, //  Br_H
	2.2, //   I_H
	1.2, // Met_D
};

//! Returns true if the XScore atom type is hydrophobic.
inline bool is_hydrophobic(const size_t t)
{
	return t ==  0 || t == 10 || t == 11 || t == 12 || t == 13;
}

//! Returns true if the XScore atom type is a hydrogen bond donor.
inline bool is_hbdonor(const size_t t)
{
	return t ==  3 || t ==  5 || t ==  7 || t == 14;
}

//! Returns true if the XScore atom type is a hydrogen bond acceptor.
inline bool is_hbacceptor(const size_t t)
{
	return t ==  4 || t ==  5 || t ==  6 || t ==  7;
}

//! Returns true if the two XScore atom types are a pair of hydrogen bond donor and acceptor.
inline bool is_hbond(const size_t t0, const size_t t1)
{
	return (is_hbdonor(t0) && is_hbacceptor(t1)) || (is_hbdonor(t1) && is_hbacceptor(t0));
}

double scoring_function::score(const size_t t0, const size_t t1, const double r)
{
	assert(r <= Cutoff);

	// Calculate the surface distance d.
	const double d = r - (vdw[t0] + vdw[t1]);

	// The scoring function is a weighted sum of 5 terms.
	// The first 3 terms depend on d only, while the latter 2 terms depend on t0, t1 and d.
	return (-0.035579) * exp(-4 * d * d)
		+  (-0.005156) * exp(-0.25 * (d - 3.0) * (d - 3.0))
		+  ( 0.840245) * (d > 0 ? 0.0 : d * d)
		+  (-0.035069) * ((is_hydrophobic(t0) && is_hydrophobic(t1)) ? ((d >= 1.5) ? 0.0 : ((d <= 0.5) ? 1.0 : 1.5 - d)) : 0.0)
		+  (-0.587439) * ((is_hbond(t0, t1)) ? ((d >= 0) ? 0.0 : ((d <= -0.7) ? 1 : d * (-1.4285714285714286))): 0.0);
}

void scoring_function::score5(float* const v, const size_t t0, const size_t t1, const double r2)
{
	assert(r2 <= Cutoff_Sqr);

	// Calculate the surface distance d.
	const double d = sqrt(r2) - (vdw[t0] + vdw[t1]);

	// The scoring function is a weighted sum of 5 terms.
	// The first 3 terms depend on d only, while the latter 2 terms depend on t0, t1 and d.
	v[0] += exp(-4 * d * d);
	v[1] += exp(-0.25 * (d - 3.0) * (d - 3.0));
	v[2] += (d > 0 ? 0.0 : d * d);
	v[3] += ((is_hydrophobic(t0) && is_hydrophobic(t1)) ? ((d >= 1.5) ? 0.0 : ((d <= 0.5) ? 1.0 : 1.5 - d)) : 0.0);
	v[4] += ((is_hbond(t0, t1)) ? ((d >= 0) ? 0.0 : ((d <= -0.7) ? 1 : d * (-1.4285714285714286))): 0.0);
}

void scoring_function::precalculate(const size_t t0, const size_t t1, const vector<double>& rs)
{
	vector<scoring_function_element>& p = (*this)[triangular_matrix_restrictive_index(t0, t1)];
	assert(p.size() == Num_Samples);

	// Calculate the value of scoring function evaluated at (t0, t1, d).
	for (size_t i = 0; i < Num_Samples; ++i)
	{
		p[i].e = score(t0, t1, rs[i]);
	}

	// Calculate the dor of scoring function evaluated at (t0, t1, d).
	for (size_t i = 1; i < Num_Samples - 1; ++i)
	{
		p[i].dor = (p[i + 1].e - p[i].e) / ((rs[i + 1] - rs[i]) * rs[i]);
	}
	p.front().dor = 0;
	p.back().dor = 0;
}

scoring_function_element scoring_function::evaluate(const size_t type_pair_index, const double r2) const
{
	assert(r2 <= Cutoff_Sqr);
	return (*this)[type_pair_index][static_cast<size_t>(Factor * r2)];
}
