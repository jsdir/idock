#pragma once
#ifndef IDOCK_SOURCE_CL_HPP
#define IDOCK_SOURCE_CL_HPP

#include <vector>
using namespace std;

/// Hardcodes idock.cl into a char vector.
class source : public vector<char>
{
public:
	source();
};

#endif
