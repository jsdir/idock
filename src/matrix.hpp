#pragma once
#ifndef IDOCK_MATRIX_HPP
#define IDOCK_MATRIX_HPP

//! Returns the flattened 1D index of a 2D index (i, j) where j is the lowest dimension.
inline size_t mr(const size_t i, const size_t j)
{
	assert(i <= j);
	return i + j * (j + 1) / 2;
}

//! Returns the flattened 1D index of a 2D index (i, j) where either i or j is the lowest dimension.
inline size_t mp(const size_t i, const size_t j)
{
	return (i <= j) ? mr(i, j) : mr(j, i);
}

#endif
