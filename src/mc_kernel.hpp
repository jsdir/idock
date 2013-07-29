#pragma once
#ifndef IDOCK_MC_KERNEL_HPP
#define IDOCK_MC_KERNEL_HPP

#include <vector>
using namespace std;

static const size_t sf_n = 15;

class mc_kernel
{
public:
	// Update grid map pointers.
	virtual void update(const vector<vector<float> > h_maps, const vector<size_t>& xs) = 0;

	// Launch the kernel.
	virtual void launch(vector<float>& h_ex, const vector<int>& h_lig, const int nv, const int nf, const int na, const int np) = 0;

	// Free memory objects.
	virtual ~mc_kernel() {};
};

#endif
