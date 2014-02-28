#pragma once
#ifndef IDOCK_ARRAY_HPP
#define IDOCK_ARRAY_HPP

#include <array>
#include "vec3.hpp"

const std::array<double, 4> qtn4id{1, 0, 0, 0}; ///< Identity quaternion.

/// Transforms a vector by a 3x3 matrix.
vec3 operator*(const std::array<double, 9>& m, const vec3& v);

std::array<double, 4> vec4_to_qtn4(const vec3& axis, const double angle);
std::array<double, 4> vec3_to_qtn4(const vec3& rotation);
double norm_sqr(const std::array<double, 4>& a);
double norm(const std::array<double, 4>& a);
bool normalized(const std::array<double, 4>& a);
std::array<double, 4> normalize(const std::array<double, 4>& a);
std::array<double, 9> qtn4_to_mat3(const std::array<double, 4>& a);

/// Returns the product of two quaternions.
std::array<double, 4> operator*(const std::array<double, 4>& q1, const std::array<double, 4>& q2);

#endif
