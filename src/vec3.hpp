#pragma once
#ifndef IDOCK_VEC3_HPP
#define IDOCK_VEC3_HPP

#include <vector>
#include <boost/array.hpp>
using std::vector;
using boost::array;

/// Returns the square of a generic value.
template<typename T>
inline T sqr(const T x)
{
	return x * x;
}

/// Returns true if the absolute difference between two floating point values is within the constant tolerance.
inline bool eq(const double a, const double b)
{
	return fabs(a - b) < 1e-5;
}

/// Represents a vector of 3 floating point elements.
class vec3 : private array<double, 3>
{
public:
	/// Constructs a vector with uninitialized values.
	vec3() {}

	/// Constructs a vector with specified values.
	vec3(const double d0, const double d1, const double d2)
	{
		elems[0] = d0;
		elems[1] = d1;
		elems[2] = d2;
	}

	/// Assigns a value to all the 3 elements.
	void assign(const double s)
	{
		elems[0] = s;
		elems[1] = s;
		elems[2] = s;
	}

	/// Returns a constant reference to the element at index i.
	const double& operator[](const size_t i) const
	{
		BOOST_ASSERT(i < 3);
		return elems[i];
	}

	/// Returns a mutable reference to the element at index i.
	double& operator[](const size_t i)
	{
		BOOST_ASSERT(i < 3);
		return elems[i];
	}

	/// Returns true is the vector is (0, 0, 0).
	bool zero() const
	{
		return (eq(elems[0], 0) && eq(elems[1], 0) && eq(elems[2], 0));
	}

	/// Returns the square norm.
	double norm_sqr() const
	{
		return sqr(elems[0]) + sqr(elems[1]) + sqr(elems[2]);
	}

	/// Returns the norm.
	double norm() const
	{
		return sqrt(norm_sqr());
	}

	/// Returns true if the norm equals 1.
	bool normalized() const
	{
		return eq(norm_sqr(), 1);
	}

	/// Normalize the vector.
	vec3 normalize() const
	{
		const double factor = 1 / norm();
		return vec3(factor * elems[0], factor * elems[1], factor * elems[2]);
	}

	/// Returns the dot product of the current vector and the given vector.
	double operator*(const vec3& v) const
	{
		return elems[0] * v[0] + elems[1] * v[1] + elems[2] * v[2];
	}

	/// Returns the result of pairwise multiplication of the current vector and the given vector.
	vec3 operator*(const array<size_t, 3>& v) const
	{
		return vec3(elems[0] * v[0], elems[1] * v[1], elems[2] * v[2]);
	}

	/// Returns the result of pairwise addition of the current vector and the given vector.
	vec3 operator+(const vec3& v) const
	{
		return vec3(elems[0] + v[0], elems[1] + v[1], elems[2] + v[2]);
	}

	/// Returns the result of pairwise subtraction of the current vector and the given vector.
	vec3 operator-(const vec3& v) const
	{
		return vec3(elems[0] - v[0], elems[1] - v[1], elems[2] - v[2]);
	}

	/// Pairwise add a given vector to the current vector.
	const vec3& operator+=(const vec3& v)
	{
		elems[0] += v[0];
		elems[1] += v[1];
		elems[2] += v[2];
		return *this;
	}

	/// Pairwise subtract a given vector from the current vector.
	const vec3& operator-=(const vec3& v)
	{
		elems[0] -= v[0];
		elems[1] -= v[1];
		elems[2] -= v[2];
		return *this;
	}
};

const vec3 zero3(0, 0, 0); ///< Constant vector with all the 3 elements of zero.

/// Pairwise multiply a constant to the current vector.
inline vec3 operator*(const double s, const vec3& v)
{
	return vec3(s * v[0], s * v[1], s * v[2]);
}

/// Returns the normalized vector of a vector.
inline vec3 normalize(const vec3& v)
{
	return (1 / v.norm()) * v;
}

/// Returns the cross product of two vectors.
inline vec3 cross_product(const vec3& a, const vec3& b)
{
	return vec3(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

/// Returns the square distance between two vectors.
inline double distance_sqr(const vec3& a, const vec3& b)
{
	return sqr(a[0] - b[0]) + sqr(a[1] - b[1]) + sqr(a[2] - b[2]);
}

/// Returns the accumulated square distance between two vectors of vectors.
inline double distance_sqr(const vector<vec3>& a, const vector<vec3>& b)
{
	const size_t n = a.size();
	BOOST_ASSERT(n > 0);
	BOOST_ASSERT(n == b.size());
	double sum = 0;
	for (size_t i = 0; i < n; ++i)
		sum += distance_sqr(a[i], b[i]);
	return sum;
}

#endif
