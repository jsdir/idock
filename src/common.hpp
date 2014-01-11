#pragma once
#ifndef IDOCK_COMMON_HPP
#define IDOCK_COMMON_HPP

#include <vector>
#include <string>

// These classes are widely used across the entire program.
using std::vector;
using std::string;

/// igrow uses double precision floating point computation by default.
/// This could possible be demoted to single precision for better performance.
typedef double fl;

const fl epsilon = static_cast<fl>(0.00001); ///< Tolerance for equality comparison of two floating point values.

/// Returns true if the absolute difference between two floating point values is within the constant tolerance.
inline bool eq(const fl a, const fl b)
{
	return fabs(a - b) < epsilon;
}

#endif
