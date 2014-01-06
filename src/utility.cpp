#include <cmath>
#include <cassert>
#include "utility.hpp"

//! Returns the flattened 1D index of a triangular 2D index (x, y) where x is the lowest dimension.
size_t mr(const size_t x, const size_t y)
{
	assert(x <= y);
	return (y*(y+1)>>1) + x;
}

//! Returns the flattened 1D index of a triangular 2D index (x, y) where either x or y is the lowest dimension.
size_t mp(const size_t x, const size_t y)
{
	return x <= y ? mr(x, y) : mr(y, x);
}

//! Returns the square norm.
float norm_sqr(const array<float, 3>& a)
{
	return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

//! Returns the norm.
float norm(const array<float, 3>& a)
{
	return sqrt(norm_sqr(a));
}

//! Returns true if the norm equals 1.
bool normalized(const array<float, 3>& a)
{
	return norm_sqr(a) - 1.0f < 3e-3f;
}

//! Normalize the vector.
array<float, 3> normalize(const array<float, 3>& a)
{
	const float norm_inv = 1.0f / norm(a);
	return
	{
		a[0] * norm_inv,
		a[1] * norm_inv,
		a[2] * norm_inv,
	};
}

//! Returns pairwise addition of 2 given arrays.
array<float, 3> operator+(const array<float, 3>& a, const array<float, 3>& b)
{
	return
	{
		a[0] + b[0],
		a[1] + b[1],
		a[2] + b[2],
	};
}

//! Returns pairwise subtraction of 2 given arrays.
array<float, 3> operator-(const array<float, 3>& a, const array<float, 3>& b)
{
	return
	{
		a[0] - b[0],
		a[1] - b[1],
		a[2] - b[2],
	};
}

//! Pairwise add a given vector to the current vector.
void operator+=(array<float, 3>& a, const array<float, 3>& b)
{
	a[0] += b[0];
	a[1] += b[1];
	a[2] += b[2];
}

//! Pairwise subtract a given vector from the current vector.
void operator-=(array<float, 3>& a, const array<float, 3>& b)
{
	a[0] -= b[0];
	a[1] -= b[1];
	a[2] -= b[2];
}

//! Pairwise multiply a constant to an array.
array<float, 3> operator*(const float s, const array<float, 3>& a)
{
	return
	{
		s * a[0],
		s * a[1],
		s * a[2],
	};
}

//! Returns the cross product of two vectors.
array<float, 3> operator*(const array<float, 3>& a, const array<float, 3>& b)
{
	return
	{
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0],
	};
}

//! Returns the square distance between two arrays.
float distance_sqr(const array<float, 3>& a, const array<float, 3>& b)
{
	const float d0 = a[0] - b[0];
	const float d1 = a[1] - b[1];
	const float d2 = a[2] - b[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

//! Returns the square norm of current quaternion.
float norm_sqr(const array<float, 4>& a)
{
	return a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3];
}

//! Returns the norm of current quaternion.
float norm(const array<float, 4>& a)
{
	return sqrt(norm_sqr(a));
}

//! Returns true if the current quaternion is normalized.
bool normalized(const array<float, 4>& a)
{
	return fabs(norm_sqr(a) - 1.0f) < 3e-3f;
}

//! Returns a normalized quaternion of current quaternion.
array<float, 4> normalize(const array<float, 4>& a)
{
	const float norm_inv = 1.0f / norm(a);
	return
	{
		a[0] * norm_inv,
		a[1] * norm_inv,
		a[2] * norm_inv,
		a[3] * norm_inv,
	};
}

//! Constructs a quaternion by a normalized axis and a rotation angle.
array<float, 4> vec4_to_qtn4(const array<float, 3>& axis, const float angle)
{
	assert(normalized(axis));
	const float h = angle * 0.5f;
	const float s = sin(h);
	const float c = cos(h);
	return
	{
		c,
		s * axis[0],
		s * axis[1],
		s * axis[2],
	};
}

//! Returns the product of two quaternions.
array<float, 4> operator*(const array<float, 4>& a, const array<float, 4>& b)
{
	assert(normalized(a));
	assert(normalized(b));
	return
	{
		a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
		a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
		a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
		a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
	};
}

//! Transforms the current quaternion into a 3x3 transformation matrix, e.g. quaternion(1, 0, 0, 0) => identity matrix.
array<float, 9> qtn4_to_mat3(const array<float, 4>& a)
{
	assert(normalized(a));
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
	return
	{
		ww+xx-yy-zz, 2*(-wz+xy), 2*(wy+xz),
		2*(wz+xy), ww-xx+yy-zz, 2*(-wx+yz),
		2*(-wy+xz), 2*(wx+yz), ww-xx-yy+zz,
	};
}

//! Transforms a vector by a 3x3 matrix.
array<float, 3> operator*(const array<float, 9>& m, const array<float, 3>& v)
{
	return
	{
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
		m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
		m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
	};
}

void safe_function::operator()(function<void(void)>&& f)
{
	lock_guard<mutex> guard(m);
	f();
}

template <typename T>
void safe_counter<T>::init(const T z)
{
	n = z;
	i = 0;
}

template <typename T>
void safe_counter<T>::increment()
{
	lock_guard<mutex> guard(m);
	if (++i == n) cv.notify_one();
}

template <typename T>
void safe_counter<T>::wait()
{
	unique_lock<mutex> lock(m);
	if (i < n) cv.wait(lock);
}

template class safe_counter<size_t>;

template <typename T>
void safe_vector<T>::safe_push_back(const T x)
{
	lock_guard<mutex> guard(m);
	this->push_back(x);
	cv.notify_one();
}

template <typename T>
T safe_vector<T>::safe_pop_back()
{
	unique_lock<mutex> lock(m);
	if (this->empty()) cv.wait(lock);
	const T x = this->back();
	this->pop_back();
	return x;
}

template class safe_vector<int>;
