#pragma once
#ifndef IDOCK_UTILITY_HPP
#define IDOCK_UTILITY_HPP

#include <array>
#include <vector>
#include <cmath>
using namespace std;

/// Returns the flattened 1D index of a triangular 2D index (x, y) where x is the lowest dimension.
inline size_t mr(const size_t x, const size_t y)
{
	return (y*(y+1)>>1) + x;
}

/// Returns the flattened 1D index of a triangular 2D index (x, y) where either x or y is the lowest dimension.
inline size_t mp(const size_t x, const size_t y)
{
	return x <= y ? mr(x, y) : mr(y, x);
}

/// Returns an array containing 3 given floats.
inline array<float, 3> make_array(const float d0, const float d1, const float d2)
{
	array<float, 3> a;
	a[0] = d0;
	a[1] = d1;
	a[2] = d2;
	return a;
}

const array<float, 3> zero3 = make_array(0.0f, 0.0f, 0.0f); ///< Constant vector with all the 3 elements of zero.

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

/// Returns pairwise addition of 2 given arrays.
inline array<float, 3> operator+(const array<float, 3>& a, const array<float, 3>& b)
{
	return make_array
	(
		a[0] + b[0],
		a[1] + b[1],
		a[2] + b[2]
	);
}

/// Returns pairwise subtraction of 2 given arrays.
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

/// Pairwise multiply a constant to an array.
inline array<float, 3> operator*(const float s, const array<float, 3>& a)
{
	return make_array
	(
		s * a[0],
		s * a[1],
		s * a[2]
	);
}

/// Returns pairwise multiplication of 2 given arrays.
inline array<float, 3> operator*(const array<float, 3>& a, const array<size_t, 3>& b)
{
	return make_array
	(
		a[0] * b[0],
		a[1] * b[1],
		a[2] * b[2]
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

/// Returns the square distance between two arrays.
inline float distance_sqr(const array<float, 3>& a, const array<float, 3>& b)
{
	const float d0 = a[0] - b[0];
	const float d1 = a[1] - b[1];
	const float d2 = a[2] - b[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

/// Returns an array containing 4 given floats.
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
//	assert(normalized(axis));
	const float h = angle * 0.5f;
	const float s = sin(h);
	const float c = cos(h);
	return make_array
	(
		c,
		s * axis[0],
		s * axis[1],
		s * axis[2]
	);
}

/// Constructs a quaternion by a rotation vector.
inline array<float, 4> vec3_to_qtn4(const array<float, 3>& rotation)
{
	if (norm_sqr(rotation) < 1e-5f)
	{
		return make_array(1.0f, 0.0f, 0.0f, 0.0f);
	}
	else
	{
		const float h = norm(rotation) * 0.5;
		const array<float, 3> axis = (0.5f / h) * rotation;
//		assert(normalized(axis));
		const float s = sin(h);
		const float c = cos(h);
		return make_array
		(
			c,
			s * axis[0],
			s * axis[1],
			s * axis[2]
		);
	}
}

/// Returns the product of two quaternions.
inline array<float, 4> operator*(const array<float, 4>& a, const array<float, 4>& b)
{
//	assert(normalized(a));
//	assert(normalized(b));
    return make_array
	(
		a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
		a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
		a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
		a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]
	);
}

/// Returns an array containing 9 given floats.
inline array<float, 9> make_array(const float d0, const float d1, const float d2, const float d3, const float d4, const float d5, const float d6, const float d7, const float d8)
{
	array<float, 9> a;
	a[0] = d0;
	a[1] = d1;
	a[2] = d2;
	a[3] = d3;
	a[4] = d4;
	a[5] = d5;
	a[6] = d6;
	a[7] = d7;
	a[8] = d8;
	return a;
}

/// Transforms the current quaternion into a 3x3 transformation matrix, e.g. quaternion(1, 0, 0, 0) => identity matrix.
inline array<float, 9> qtn4_to_mat3(const array<float, 4>& a)
{
//	assert(normalized(a));
	const float ww = a[0]*a[0];
	const float wx = a[0]*a[1];
	const float wy = a[0]*a[2];
	const float wz = a[0]*a[3];
	const float xx = a[1]*a[1];
	const float xy = a[1]*a[2];
	const float xz = a[1]*a[3];
	const float yy = a[2]*a[2];
	const float yz = a[2]*a[3];
	const float zz = a[3]*a[3];

	// http://www.boost.org/doc/libs/1_46_1/libs/math/quaternion/TQE.pdf
	// http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	return make_array
	(
		ww+xx-yy-zz, 2*(-wz+xy), 2*(wy+xz),
		2*(wz+xy), ww-xx+yy-zz, 2*(-wx+yz),
		2*(-wy+xz), 2*(wx+yz), ww-xx-yy+zz
	);
}

/// Transforms a vector by a 3x3 matrix.
inline array<float, 3> operator*(const array<float, 9>& m, const array<float, 3>& v)
{
	return make_array
	(
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
		m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
		m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
	);
}

#endif
