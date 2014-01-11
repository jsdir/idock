#pragma once
#ifndef IDOCK_UTILITY_HPP
#define IDOCK_UTILITY_HPP

#include <vector>
#include <functional>
#include <condition_variable>
using namespace std;

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
