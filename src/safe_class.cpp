#include "safe_class.hpp"

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
