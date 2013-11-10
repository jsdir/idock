#pragma once
#ifndef IDOCK_CU_MC_KERNEL_HPP
#define IDOCK_CU_MC_KERNEL_HPP

static const size_t sf_n = 15;

class cu_mc_kernel
{
public:
	cu_mc_kernel(const int device_id) : device_id(device_id)
	{
	}

	// Initialize kernel.
	int initialize(const int tid, const vector<float>& h_sf_e, const vector<float>& h_sf_d, const size_t h_sf_ns, const float* h_corner0, const float* h_corner1, const int* h_num_probes, const float h_granularity_inverse, const int num_mc_tasks, const int h_ng, const unsigned long h_seed);

	// Update grid map pointers.
	int update(const int tid, const vector<vector<float> > h_maps, const vector<size_t>& xs);

	// Launch the kernel.
	void launch(vector<float>& h_ex, const vector<int>& h_lig, const int nv, const int nf, const int na, const int np);

	// Free device memory.
	~cu_mc_kernel();
private:
	const int device_id;
	float* d_sf_e;
	float* d_sf_d;
	float* d_maps[sf_n];
	int num_mc_tasks;
};

#endif
