#pragma once
#ifndef IDOCK_KERNEL_HPP
#define IDOCK_KERNEL_HPP

#include <array>
using namespace std;

void monte_carlo(float* const s0e, const int* const lig, const int nv, const int nf, const int na, const int np, const int nbi, const float* const sfe, const float* const sfd, const int sfs, const array<float, 3> cr0, const array<float, 3> cr1, const array<int, 3> npr, const float gri, const vector<vector<float>>& mps, const int gid, const int gds);

#endif
