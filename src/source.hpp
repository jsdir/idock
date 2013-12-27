#pragma once
#ifndef IDOCK_SOURCE_HPP
#define IDOCK_SOURCE_HPP

#include <vector>
using namespace std;

//! Represents the kernel source as a char vector.
class source : public vector<char>
{
public:
	explicit source();
};

#endif
