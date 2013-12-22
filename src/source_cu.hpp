#pragma once
#ifndef IDOCK_SOURCE_CU_HPP
#define IDOCK_SOURCE_CU_HPP

#include <vector>
using namespace std;

/// Hardcodes idock.fatbin into a char vector.
class source : public vector<unsigned char>
{
public:
	source();
};

#endif
