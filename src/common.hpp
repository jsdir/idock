#pragma once
#ifndef IDOCK_COMMON_HPP
#define IDOCK_COMMON_HPP

#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/random.hpp>
#include <boost/assert.hpp>

// These classes are widely used across the entire program.
using std::vector;
using std::string;

/// igrow uses double precision floating point computation by default.
/// This could possible be demoted to single precision for better performance.
typedef double fl;

// Choose the appropriate Mersenne Twister engine for random number generation on 32-bit or 64-bit platform.
#if defined(__x86_64) || defined(__x86_64__) || defined(__amd64) || defined(__amd64__) || defined(_M_X64) || defined(_M_AMD64)
typedef boost::random::mt19937_64 mt19937eng;
#else
typedef boost::random::mt19937 mt19937eng;
#endif

const fl epsilon = static_cast<fl>(0.00001); ///< Tolerance for equality comparison of two floating point values.

/// Returns true if the absolute difference between two floating point values is within the constant tolerance.
inline bool eq(const fl a, const fl b)
{
	return fabs(a - b) < epsilon;
}

#endif
