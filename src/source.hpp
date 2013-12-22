#pragma once
#ifndef IDOCK_SOURCE_HPP
#define IDOCK_SOURCE_HPP

#include <vector>
using namespace std;

/// Represents the kernel source as a char vector.
template <typename T>
class source : public vector<T>
{
public:
	source();
};

#endif
