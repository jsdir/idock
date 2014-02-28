#include <cmath>
#include <boost/assert.hpp>
#include "array.hpp"

inline bool zero(const double a)
{
	return fabs(a) < 1e-5;
}

/// Returns true if the absolute difference between two floating point values is within the constant tolerance.
inline bool eq(const double a, const double b)
{
	return zero(a - b);
}

/// Returns true is the vector is (0, 0, 0).
bool zero(const array<double, 3>& a)
{
	return zero(a[0]) && zero(a[1]) && zero(a[2]);
}

/// Returns the square norm.
double norm_sqr(const array<double, 3>& a)
{
	return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
}

/// Returns the norm.
double norm(const array<double, 3>& a)
{
	return sqrt(norm_sqr(a));
}

/// Returns true if the norm equals 1.
bool normalized(const array<double, 3>& a)
{
	return eq(norm_sqr(a), 1);
}

/// Returns the dot product of the current vector and the given vector.
double operator*(const array<double, 3>& t, const array<double, 3>& v)
{
	return t[0] * v[0] + t[1] * v[1] + t[2] * v[2];
}

/// Returns the result of pairwise multiplication of the current vector and the given vector.
array<double, 3> operator*(const array<double, 3>& t, const array<size_t, 3>& v)
{
	return
	{
		t[0] * v[0],
		t[1] * v[1],
		t[2] * v[2],
	};
}

/// Returns the result of pairwise addition of the current vector and the given vector.
array<double, 3> operator+(const array<double, 3>& t, const array<double, 3>& v)
{
	return
	{
		t[0] + v[0],
		t[1] + v[1],
		t[2] + v[2],
	};
}

/// Returns the result of pairwise subtraction of the current vector and the given vector.
array<double, 3> operator-(const array<double, 3>& t, const array<double, 3>& v)
{
	return
	{
		t[0] - v[0],
		t[1] - v[1],
		t[2] - v[2],
	};
}

/// Pairwise add a given vector to the current vector.
void operator+=(array<double, 3>& t, const array<double, 3>& v)
{
	t[0] += v[0];
	t[1] += v[1];
	t[2] += v[2];
}

/// Pairwise subtract a given vector from the current vector.
void operator-=(array<double, 3>& t, const array<double, 3>& v)
{
	t[0] -= v[0];
	t[1] -= v[1];
	t[2] -= v[2];
}

/// Pairwise multiply a constant to the current vector.
array<double, 3> operator*(const double s, const array<double, 3>& v)
{
	return
	{
		s * v[0],
		s * v[1],
		s * v[2],
	};
}

/// Returns the normalized vector of a vector.
array<double, 3> normalize(const array<double, 3>& v)
{
	return (1 / norm(v)) * v;
}

/// Returns the cross product of two vectors.
array<double, 3> cross_product(const array<double, 3>& a, const array<double, 3>& b)
{
	return
	{
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0],
	};
}

/// Returns the square distance between two vectors.
double distance_sqr(const array<double, 3>& a, const array<double, 3>& b)
{
	const auto d0 = a[0] - b[0];
	const auto d1 = a[1] - b[1];
	const auto d2 = a[2] - b[2];
	return d0*d0 + d1*d1 + d2*d2;
}

/// Returns the accumulated square distance between two vectors of vectors.
double distance_sqr(const vector<array<double, 3>>& a, const vector<array<double, 3>>& b)
{
	const size_t n = a.size();
	BOOST_ASSERT(n > 0);
	BOOST_ASSERT(n == b.size());
	double sum = 0;
	for (size_t i = 0; i < n; ++i)
	{
		sum += distance_sqr(a[i], b[i]);
	}
	return sum;
}

/// Transforms a vector by a 3x3 matrix.
array<double, 3> operator*(const array<double, 9>& m, const array<double, 3>& v)
{
	return
	{
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
		m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
		m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
	};
}

array<double, 4> vec4_to_qtn4(const array<double, 3>& axis, const double angle)
{
	BOOST_ASSERT(normalized(axis));
	const double h = angle * 0.5;
	const double s = sin(h);
	const double c = cos(h);
	return
	{
		c,
		s * axis[0],
		s * axis[1],
		s * axis[2],
	};
}

array<double, 4> vec3_to_qtn4(const array<double, 3>& rotation)
{
	if (zero(rotation))
	{
		return qtn4id;
	}
	else
	{
		const auto angle = norm(rotation);
		const auto axis = (1 / angle) * rotation;
		return vec4_to_qtn4(axis, angle);
	}
}

double norm_sqr(const array<double, 4>& a)
{
	return a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3];
}

double norm(const array<double, 4>& a)
{
	return sqrt(norm_sqr(a));
}

bool normalized(const array<double, 4>& a)
{
	return fabs(norm_sqr(a) - 1.0) < 1e-5;
}

array<double, 4> normalize(const array<double, 4>& a)
{
	const auto norm_inv = 1 / norm(a);
	return
	{
		norm_inv * a[0],
		norm_inv * a[1],
		norm_inv * a[2],
		norm_inv * a[3],
	};
}

array<double, 9> qtn4_to_mat3(const array<double, 4>& a)
{
	assert(normalized(a));
	const auto ww = a[0]*a[0];
	const auto wx = a[0]*a[1];
	const auto wy = a[0]*a[2];
	const auto wz = a[0]*a[3];
	const auto xx = a[1]*a[1];
	const auto xy = a[1]*a[2];
	const auto xz = a[1]*a[3];
	const auto yy = a[2]*a[2];
	const auto yz = a[2]*a[3];
	const auto zz = a[3]*a[3];

	// http://www.boost.org/doc/libs/1_46_1/libs/math/quaternion/TQE.pdf
	// http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	return
	{
		ww+xx-yy-zz, 2*(-wz+xy), 2*(wy+xz),
		2*(wz+xy), ww-xx+yy-zz, 2*(-wx+yz),
		2*(-wy+xz), 2*(wx+yz), ww-xx-yy+zz,
	};
}

/// Returns the product of two quaternions.
array<double, 4> operator*(const array<double, 4>& q1, const array<double, 4>& q2)
{
    return
	{
		q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3],
		q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2],
		q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1],
		q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0],
	};
}
