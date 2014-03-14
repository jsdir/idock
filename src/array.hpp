#pragma once
#ifndef IDOCK_ARRAY_HPP
#define IDOCK_ARRAY_HPP

#include <array>
#include <vector>
using namespace std;

const array<double, 4> qtn4id{1, 0, 0, 0}; //!< Identity quaternion.

double norm_sqr(const array<double, 3>& a);
double norm(const array<double, 3>& a);
bool normalized(const array<double, 3>& a);
array<double, 3> normalize(const array<double, 3>& v);
double operator*(const array<double, 3>& t, const array<double, 3>& v);
array<double, 3> operator*(const array<double, 3>& t, const array<size_t, 3>& v);
array<double, 3> operator+(const array<double, 3>& t, const array<double, 3>& v);
array<double, 3> operator-(const array<double, 3>& t, const array<double, 3>& v);
void operator+=(array<double, 3>& t, const array<double, 3>& v);
void operator-=(array<double, 3>& t, const array<double, 3>& v);
array<double, 3> operator*(const double s, const array<double, 3>& v);
array<double, 3> cross_product(const array<double, 3>& a, const array<double, 3>& b);
double distance_sqr(const array<double, 3>& a, const array<double, 3>& b);
double distance_sqr(const vector<array<double, 3>>& a, const vector<array<double, 3>>& b);

//! Transforms a vector by a 3x3 matrix.
array<double, 3> operator*(const array<double, 9>& m, const array<double, 3>& v);

array<double, 4> vec4_to_qtn4(const array<double, 3>& axis, const double angle);
array<double, 4> vec3_to_qtn4(const array<double, 3>& rotation);
double norm_sqr(const array<double, 4>& a);
double norm(const array<double, 4>& a);
bool normalized(const array<double, 4>& a);
array<double, 4> normalize(const array<double, 4>& a);
array<double, 9> qtn4_to_mat3(const array<double, 4>& a);

//! Returns the product of two quaternions.
array<double, 4> operator*(const array<double, 4>& a, const array<double, 4>& b);

#endif
