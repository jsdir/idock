#pragma once
#ifndef IDOCK_KERNEL_HPP
#define IDOCK_KERNEL_HPP

const size_t sf_n = 15;

class kernel
{
public:
	// Initialize kernel.
	explicit kernel(const float* h_sf_e, const float* h_sf_d, const int h_sf_ns, const int h_sf_ne, const float* h_corner0, const float* h_corner1, const float* h_num_probes, const float h_granularity_inverse, const int num_mc_tasks, const int num_generations);

	// Update grid map pointers.
	void update(const float* h_maps[], const size_t num_probes_product, const size_t* xs, const size_t n);

	// Launch the kernel.
	void launch(float* h_s0, const float* lig, const int n, const size_t oz, const size_t og, const size_t* seed);

	// Free device memory.
	~kernel();
private:
	float* d_sf_e;
	float* d_sf_d;
	float* d_maps[sf_n];
	int num_mc_tasks;
	int num_generations;
};

#endif
