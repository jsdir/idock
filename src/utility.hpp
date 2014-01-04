#pragma once
#ifndef IDOCK_UTILITY_HPP
#define IDOCK_UTILITY_HPP

#include <array>
#include <vector>
#include <functional>
#include <condition_variable>
using namespace std;

//! Returns the flattened 1D index of a triangular 2D index (x, y) where x is the lowest dimension.
size_t mr(const size_t x, const size_t y);

//! Returns the flattened 1D index of a triangular 2D index (x, y) where either x or y is the lowest dimension.
size_t mp(const size_t x, const size_t y);

//! Returns the square norm.
float norm_sqr(const array<float, 3>& a);

//! Returns the norm.
float norm(const array<float, 3>& a);

//! Returns true if the norm equals 1.
bool normalized(const array<float, 3>& a);

//! Normalize the vector.
array<float, 3> normalize(const array<float, 3>& a);

//! Returns pairwise addition of 2 given arrays.
array<float, 3> operator+(const array<float, 3>& a, const array<float, 3>& b);

//! Returns pairwise subtraction of 2 given arrays.
array<float, 3> operator-(const array<float, 3>& a, const array<float, 3>& b);

//! Pairwise add a given vector to the current vector.
void operator+=(array<float, 3>& a, const array<float, 3>& b);

//! Pairwise subtract a given vector from the current vector.
void operator-=(array<float, 3>& a, const array<float, 3>& b);

//! Pairwise multiply a constant to an array.
array<float, 3> operator*(const float s, const array<float, 3>& a);

//! Returns the cross product of two vectors.
array<float, 3> operator*(const array<float, 3>& a, const array<float, 3>& b);

//! Returns the square distance between two arrays.
float distance_sqr(const array<float, 3>& a, const array<float, 3>& b);

//! Returns the square norm of current quaternion.
float norm_sqr(const array<float, 4>& a);

//! Returns the norm of current quaternion.
float norm(const array<float, 4>& a);

//! Returns true if the current quaternion is normalized.
bool normalized(const array<float, 4>& a);

//! Returns a normalized quaternion of current quaternion.
array<float, 4> normalize(const array<float, 4>& a);

//! Constructs a quaternion by a normalized axis and a rotation angle.
array<float, 4> vec4_to_qtn4(const array<float, 3>& axis, const float angle);

//! Returns the product of two quaternions.
array<float, 4> operator*(const array<float, 4>& a, const array<float, 4>& b);

//! Transforms the current quaternion into a 3x3 transformation matrix, e.g. quaternion(1, 0, 0, 0) => identity matrix.
array<float, 9> qtn4_to_mat3(const array<float, 4>& a);

//! Transforms a vector by a 3x3 matrix.
array<float, 3> operator*(const array<float, 9>& m, const array<float, 3>& v);

//! Represents a thread safe function.
class safe_function
{
public:
	//! Executes a function in a thread safe manner.
	void operator()(function<void(void)>&& f);
private:
	mutex m;
};

//! Represents a thread safe counter.
template <typename T>
class safe_counter
{
public:
	//! Initializes the counter to 0 and its expected hit value to z.
	void init(const T z);

	//! Increments the counter by 1 in a thread safe manner, and wakes up the calling thread waiting on the internal mutex.
	void increment();

	//! Waits until the counter reaches its expected hit value.
	void wait();
private:
	mutex m;
	condition_variable cv;
	T n; //!< Expected hit value.
	T i; //!< Counter value.
};

//! Represents a thread safe vector.
template <typename T>
class safe_vector : public vector<T>
{
public:
	using vector<T>::vector;

	//! Pushes back a new element in a thread safe manner, and wakes up the calling thread waiting on new elements to come.
	void safe_push_back(const T x);

	//! Waits until the vector is not empty, and pops back the last element in a thread safe manner.
	T safe_pop_back();
private:
	mutex m;
	condition_variable cv;
};

#endif
