#pragma once
#ifndef IDOCK_KERNEL_HPP
#define IDOCK_KERNEL_HPP

#include <vector>
using namespace std;

const size_t sf_n = 15;

class kernel
{
public:
	// Initialize kernel.
	explicit kernel(const float* h_sf_e, const float* h_sf_d, const int h_sf_ns, const int h_sf_ne, const float* h_corner0, const float* h_corner1, const int* h_num_probes, const float h_granularity_inverse, const int num_mc_tasks, const int h_ng);

	// Update grid map pointers.
	void update(const vector<vector<float> > h_maps, const size_t map_bytes, const vector<size_t>& xs);

	// Launch the kernel.
	void launch(vector<float>& h_ex, const vector<int>& h_lig, const int nv, const int nf, const int na, const int np, const unsigned long long seed);

	// Free device memory.
	~kernel();
private:
	float* d_sf_e;
	float* d_sf_d;
	float* d_maps[sf_n];
	int num_mc_tasks;
	int ng;
};

#endif
