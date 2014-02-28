#include "array.hpp"

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

std::array<double, 4> vec4_to_qtn4(const vec3& axis, const double angle)
{
	BOOST_ASSERT(axis.normalized());
	const double h = angle * 0.5;
	const double s = sin(h);
	const double c = cos(h);
	return
	{
		c,
		s * axis[0],
		s * axis[1],
		s * axis[2],
	};
}

std::array<double, 4> vec3_to_qtn4(const vec3& rotation)
{
	if (rotation.zero())
	{
		return qtn4id;
	}
	else
	{
		const double angle = rotation.norm();
		const vec3 axis = (1 / angle) * rotation;
		return vec4_to_qtn4(axis, angle);
	}
}

double norm_sqr(const std::array<double, 4>& a)
{
	return a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3];
}

double norm(const std::array<double, 4>& a)
{
	return sqrt(norm_sqr(a));
}

bool normalized(const std::array<double, 4>& a)
{
	return fabs(norm_sqr(a) - 1.0) < 1e-5;
}

std::array<double, 4> normalize(const std::array<double, 4>& a)
{
	const double norm_inv = 1.0 / norm(a);
	return
	{
		a[0] * norm_inv,
		a[1] * norm_inv,
		a[2] * norm_inv,
		a[3] * norm_inv,
	};
}

std::array<double, 9> qtn4_to_mat3(const std::array<double, 4>& a)
{
	assert(normalized(a));
	const auto ww = a[0]*a[0];
	const auto wx = a[0]*a[1];
	const auto wy = a[0]*a[2];
	const auto wz = a[0]*a[3];
	const auto xx = a[1]*a[1];
	const auto xy = a[1]*a[2];
	const auto xz = a[1]*a[3];
	const auto yy = a[2]*a[2];
	const auto yz = a[2]*a[3];
	const auto zz = a[3]*a[3];

	// http://www.boost.org/doc/libs/1_46_1/libs/math/quaternion/TQE.pdf
	// http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	return
	{
		ww+xx-yy-zz, 2*(-wz+xy), 2*(wy+xz),
		2*(wz+xy), ww-xx+yy-zz, 2*(-wx+yz),
		2*(-wy+xz), 2*(wx+yz), ww-xx-yy+zz,
	};
}

/// Returns the product of two quaternions.
std::array<double, 4> operator*(const std::array<double, 4>& q1, const std::array<double, 4>& q2)
{
    return
	{
		q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3],
		q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2],
		q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1],
		q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0],
	};
}
