#pragma once
#ifndef IDOCK_UTILITY_HPP
#define IDOCK_UTILITY_HPP

#include <array>
#include <vector>
#include <cmath>
using namespace std;

/// Returns the flattened 1D index of a 2D index (i, j) where j is the lowest dimension.
inline size_t mr(const size_t i, const size_t j)
{
	return (j*(j+1)>>1) + i;
}

/// Returns the flattened 1D index of a 2D index (i, j) where either i or j is the lowest dimension.
inline size_t mp(const size_t i, const size_t j)
{
	return i <= j ? mr(i, j) : mr(j, i);
}

/// Returns true is the vector is (0, 0, 0).
inline bool zero(const array<float, 3>& a)
{
	return a[0] < 1e-5f && a[1] < 1e-5f && a[2] < 1e-5f;
}

inline array<float, 3> make_array(const float d0, const float d1, const float d2)
{
	array<float, 3> a;
	a[0] = d0;
	a[1] = d1;
	a[2] = d2;
	return a;
}

/// Returns the result of pairwise multiplication of the current vector and the given vector.
inline array<float, 3> operator*(const array<float, 3>& a, const array<size_t, 3>& b)
{
	return make_array
	(
		a[0] * b[0],
		a[1] * b[1],
		a[2] * b[2]
	);
}

/// Returns the result of pairwise addition of the current vector and the given vector.
inline array<float, 3> operator+(const array<float, 3>& a, const array<float, 3>& b)
{
	return make_array
	(
		a[0] + b[0],
		a[1] + b[1],
		a[2] + b[2]
	);
}

/// Returns the result of pairwise subtraction of the current vector and the given vector.
inline array<float, 3> operator-(const array<float, 3>& a, const array<float, 3>& b)
{
	return make_array
	(
		a[0] - b[0],
		a[1] - b[1],
		a[2] - b[2]
	);
}

/// Pairwise add a given vector to the current vector.
inline void operator+=(array<float, 3>& a, const array<float, 3>& b)
{
	a[0] += b[0];
	a[1] += b[1];
	a[2] += b[2];
}

/// Pairwise subtract a given vector from the current vector.
inline void operator-=(array<float, 3>& a, const array<float, 3>& b)
{
	a[0] -= b[0];
	a[1] -= b[1];
	a[2] -= b[2];
}

const array<float, 3> zero3 = make_array(0.0f, 0.0f, 0.0f); ///< Constant vector with all the 3 elements of zero.

/// Pairwise multiply a constant to the current vector.
inline array<float, 3> operator*(const float s, const array<float, 3>& v)
{
	return make_array
	(
		s * v[0],
		s * v[1],
		s * v[2]
	);
}

/// Returns the cross product of two vectors.
inline array<float, 3> cross_product(const array<float, 3>& a, const array<float, 3>& b)
{
	return make_array
	(
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0]
	);
}

/// Returns the square distance between two vectors.
inline float distance_sqr(const array<float, 3>& a, const array<float, 3>& b)
{
	const float d0 = a[0] - b[0];
	const float d1 = a[1] - b[1];
	const float d2 = a[2] - b[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

/// Returns the square norm.
inline float norm_sqr(const array<float, 3>& a)
{
	return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

/// Returns the norm.
inline float norm(const array<float, 3>& a)
{
	return sqrt(norm_sqr(a));
}

/// Returns true if the norm equals 1.
inline bool normalized(const array<float, 3>& a)
{
	return norm_sqr(a) - 1.0f < 1e-5f;
}

/// Normalize the vector.
inline array<float, 3> normalize(const array<float, 3>& a)
{
	const float norm_inv = 1.0f / norm(a);
	return make_array
	(
		a[0] * norm_inv,
		a[1] * norm_inv,
		a[2] * norm_inv
	);
}

inline array<float, 4> make_array(const float d0, const float d1, const float d2, const float d3)
{
	array<float, 4> a;
	a[0] = d0;
	a[1] = d1;
	a[2] = d2;
	a[3] = d3;
	return a;
}

/// Returns the square norm of current quaternion.
inline float norm_sqr(const array<float, 4>& a)
{
	return a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3];
}

/// Returns the norm of current quaternion.
inline float norm(const array<float, 4>& a)
{
	return sqrt(norm_sqr(a));
}

/// Returns true if the current quaternion is normalized.
inline bool normalized(const array<float, 4>& a)
{
	return norm_sqr(a) - 1.0f < 1e-5f;
}

/// Returns a normalized quaternion of current quaternion.
inline array<float, 4> normalize(const array<float, 4>& a)
{
	const float norm_inv = 1.0f / norm(a);
	return make_array
	(
		a[0] * norm_inv,
		a[1] * norm_inv,
		a[2] * norm_inv,
		a[3] * norm_inv
	);
}

/// Constructs a quaternion by a normalized axis and a rotation angle.
inline array<float, 4> vec4_to_qtn4(const array<float, 3>& axis, const float angle)
{
//	assert(axis.normalized());
	const float s = sin(angle * 0.5f);
	return make_array
	(
		cos(angle * 0.5f),
		s * axis[0],
		s * axis[1],
		s * axis[2]
	);
}

/// Constructs a quaternion by a rotation vector.
inline array<float, 4> vec3_to_qtn4(const array<float, 3>& rotation)
{
	if (zero(rotation))
	{
		return make_array(1.0f, 0.0f, 0.0f, 0.0f);
	}
	else
	{
		const float angle = norm(rotation);
		const array<float, 3> axis = (1.0f / angle) * rotation;
		return vec4_to_qtn4(axis, angle);
	}
}

/// Returns the product of two quaternions.
inline array<float, 4> qtn4_mul_qtn4(const array<float, 4>& q1, const array<float, 4>& q2)
{
    return make_array
	(
		q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3],
		q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2],
		q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1],
		q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0]
	);
}

/// Transforms the current quaternion into a 3x3 transformation matrix, e.g. quaternion(1, 0, 0, 0) => identity matrix.
inline array<float, 9> qtn4_to_mat3(const array<float, 4>& q)
{
//	assert(is_normalized());
	const float aa = q[0]*q[0];
	const float ab = q[0]*q[1];
	const float ac = q[0]*q[2];
	const float ad = q[0]*q[3];
	const float bb = q[1]*q[1];
	const float bc = q[1]*q[2];
	const float bd = q[1]*q[3];
	const float cc = q[2]*q[2];
	const float cd = q[2]*q[3];
	const float dd = q[3]*q[3];

	// http://www.boost.org/doc/libs/1_46_1/libs/math/quaternion/TQE.pdf
	// http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	const array<float, 9> r =
	{
		aa+bb-cc-dd, 2*(-ad+bc), 2*(ac+bd),
		2*(ad+bc), aa-bb+cc-dd, 2*(-ab+cd),
		2*(-ac+bd), 2*(ab+cd), aa-bb-cc+dd
	};
	return r;
}

/// Transforms a vector by a 3x3 matrix.
inline array<float, 3> mat3_mul_vec3(const array<float, 9>& m, const array<float, 3>& v)
{
	return make_array
	(
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
		m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
		m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
	);
}

#endif
