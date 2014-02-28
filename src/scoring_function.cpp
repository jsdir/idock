#include "scoring_function.hpp"

const double scoring_function::Cutoff = static_cast<double>(8);
const double scoring_function::Cutoff_Sqr = Cutoff * Cutoff;
const double scoring_function::Factor = static_cast<double>(256);
const double scoring_function::Factor_Inverse = 1 / Factor;
const size_t scoring_function::Num_Samples = static_cast<size_t>(Factor * Cutoff_Sqr) + 1;

const double xs_vdw_radii[] = ///< Van der Waals radii for XScore atom types.
{
	1.9, //  0 = XS_TYPE_C_H
	1.9, //  1 = XS_TYPE_C_P
	1.8, //  2 = XS_TYPE_N_P
	1.8, //  3 = XS_TYPE_N_D
	1.8, //  4 = XS_TYPE_N_A
	1.8, //  5 = XS_TYPE_N_DA
	1.7, //  6 = XS_TYPE_O_A
	1.7, //  7 = XS_TYPE_O_DA
	2.0, //  8 = XS_TYPE_S_P
	2.1, //  9 = XS_TYPE_P_P
	1.5, // 10 = XS_TYPE_F_H
	1.8, // 11 = XS_TYPE_Cl_H
	2.0, // 12 = XS_TYPE_Br_H
	2.2, // 13 = XS_TYPE_I_H
	1.2  // 14 = XS_TYPE_Met_D
};

/// Returns Van der Waals radius from an XScore atom type.
inline double xs_vdw_radius(const size_t xs)
{
	BOOST_ASSERT(xs < XS_TYPE_SIZE);
	return xs_vdw_radii[xs];
}

/// Returns true if the XScore atom type is hydrophobic.
inline bool xs_is_hydrophobic(const size_t xs)
{
	BOOST_ASSERT(xs < XS_TYPE_SIZE);
	return xs == XS_TYPE_C_H
		|| xs == XS_TYPE_F_H
		|| xs == XS_TYPE_Cl_H
		|| xs == XS_TYPE_Br_H
		|| xs == XS_TYPE_I_H;
}

/// Returns true if the XScore atom type is a hydrogen bond donor.
inline bool xs_is_donor(const size_t xs)
{
	BOOST_ASSERT(xs < XS_TYPE_SIZE);
	return xs == XS_TYPE_N_D
		|| xs == XS_TYPE_N_DA
		|| xs == XS_TYPE_O_DA
		|| xs == XS_TYPE_Met_D;
}

/// Returns true if the XScore atom type is a hydrogen bond acceptor.
inline bool xs_is_acceptor(const size_t xs)
{
	BOOST_ASSERT(xs < XS_TYPE_SIZE);
	return xs == XS_TYPE_N_A
		|| xs == XS_TYPE_N_DA
		|| xs == XS_TYPE_O_A
		|| xs == XS_TYPE_O_DA;
}

/// Returns true if the XScore atom type is either a hydrogen bond donor or a hydrogen bond acceptor.
inline bool xs_is_donor_acceptor(const size_t xs)
{
	BOOST_ASSERT(xs < XS_TYPE_SIZE);
	return xs_is_donor(xs) || xs_is_acceptor(xs);
}

/// Returns true if the two XScore atom types are a pair of hydrogen bond donor and acceptor.
inline bool xs_hbond(const size_t xs1, const size_t xs2)
{
	return (xs_is_donor(xs1) && xs_is_acceptor(xs2))
		|| (xs_is_donor(xs2) && xs_is_acceptor(xs1));
}

double scoring_function::score(const size_t t0, const size_t t1, const double r)
{
	BOOST_ASSERT(r <= Cutoff);

	// Calculate the surface distance d.
	const double d = r - (xs_vdw_radius(t0) + xs_vdw_radius(t1));

	// The scoring function is a weighted sum of 5 terms.
	// The first 3 terms depend on d only, while the latter 2 terms depend on t0, t1 and d.
	return (-0.035579) * exp(-sqr(d * 2))
		+  (-0.005156) * exp(-sqr((d - 3.0) * 0.5))
		+  ( 0.840245) * (d > 0 ? 0.0 : d * d)
		+  (-0.035069) * ((xs_is_hydrophobic(t0) && xs_is_hydrophobic(t1)) ? ((d >= 1.5) ? 0.0 : ((d <= 0.5) ? 1.0 : 1.5 - d)) : 0.0)
		+  (-0.587439) * ((xs_hbond(t0, t1)) ? ((d >= 0) ? 0.0 : ((d <= -0.7) ? 1 : d * (-1.4285714285714286))): 0.0);
}

void scoring_function::score5(float* const v, const size_t t0, const size_t t1, const double r2)
{
	BOOST_ASSERT(r2 <= Cutoff_Sqr);

	// Calculate the surface distance d.
	const double d = sqrt(r2) - (xs_vdw_radius(t0) + xs_vdw_radius(t1));

	// The scoring function is a weighted sum of 5 terms.
	// The first 3 terms depend on d only, while the latter 2 terms depend on t0, t1 and d.
	v[0] += exp(-sqr(d * 2));
	v[1] += exp(-sqr((d - 3.0) * 0.5));
	v[2] += (d > 0 ? 0.0 : d * d);
	v[3] += ((xs_is_hydrophobic(t0) && xs_is_hydrophobic(t1)) ? ((d >= 1.5) ? 0.0 : ((d <= 0.5) ? 1.0 : 1.5 - d)) : 0.0);
	v[4] += ((xs_hbond(t0, t1)) ? ((d >= 0) ? 0.0 : ((d <= -0.7) ? 1 : d * (-1.4285714285714286))): 0.0);
}

void scoring_function::precalculate(const size_t t0, const size_t t1, const vector<double>& rs)
{
	vector<scoring_function_element>& p = (*this)[triangular_matrix_restrictive_index(t0, t1)];
	BOOST_ASSERT(p.size() == Num_Samples);

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
	BOOST_ASSERT(r2 <= Cutoff_Sqr);
	return (*this)[type_pair_index][static_cast<size_t>(Factor * r2)];
}
