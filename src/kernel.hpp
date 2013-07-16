#pragma once
#ifndef IDOCK_KERNEL_HPP
#define IDOCK_KERNEL_HPP

#include "ligand.hpp"

class kernel
{
public:
	// Initialize kernel.
	explicit kernel(const scoring_function& sf, const receptor& rec, const size_t num_mc_tasks, const size_t num_generations);

	// Update grid map pointers.
	void update(const receptor& rec, const vector<size_t>& xs);

	// Launch the kernel.
	void launch(vector<float>& h_s0, const ligand& lig);

	// Free device memory.
	~kernel();
private:
	float* d_sfe;
	float* d_sfd;
	array<float*, scoring_function::n> d_maps;
	size_t num_mc_tasks;
	size_t num_generations;
};

#endif
