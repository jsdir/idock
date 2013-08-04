#pragma once
#ifndef IDOCK_MC_KERNEL_HPP
#define IDOCK_MC_KERNEL_HPP

#include <vector>
#include "kernel.hpp"
using namespace std;

static const size_t sf_n = 15;

class mc_kernel : public kernel
{
public:
	mc_kernel(const int device_id) : kernel(device_id)
	{
	}

	// Initialize kernel.
	virtual int initialize(const int tid, const vector<float>& h_sf_e, const vector<float>& h_sf_d, const size_t h_sf_ns, const float* h_corner0, const float* h_corner1, const int* h_num_probes, const float h_granularity_inverse, const int num_mc_tasks, const int h_ng, const unsigned long h_seed) = 0;

	// Update grid map pointers.
	virtual int update(const int tid, const vector<vector<float> > h_maps, const vector<size_t>& xs) = 0;

	// Launch the kernel.
	virtual void launch(vector<float>& h_ex, const vector<int>& h_lig, const int nv, const int nf, const int na, const int np) = 0;

	// Free memory objects.
	virtual ~mc_kernel() {};
};

#endif
