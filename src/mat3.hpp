#pragma once
#ifndef IDOCK_MAT3_HPP
#define IDOCK_MAT3_HPP

#include "vec3.hpp"

/// Transforms a vector by a 3x3 matrix.
vec3 operator*(const std::array<double, 9>& m, const vec3& v)
{
	return vec3
	(
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
		m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
		m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
	);
}

#endif
