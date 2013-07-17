#include <assert.h>
#include <cuda_runtime_api.h>
#include "kernel.hpp"

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#undef  assert
#define assert(arg)
#endif

__constant__ float* c_sf_e;
__constant__ float* c_sf_d;
__constant__ int c_sf_ns;
__constant__ float3 c_corner0;
__constant__ float3 c_corner1;
__constant__ float3 c_num_probes;
__constant__ float c_granularity_inverse;
__constant__ float* c_maps[sf_n];
__constant__ int c_num_generations;

extern __shared__ float shared[];

__device__  __noinline__// __forceinline__
bool evaluate(const float* x, float* e, float* g, float* a, float* q, float* c, float* d, float* f, float* t, const float e_upper_bound)
{
	return true;
}

__global__
//__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
void bfgs(float* __restrict__ s0e, float* __restrict__ s1e, float* __restrict__ s2e, const float* lig, const int nv, const int nf, const int na, const int np, const int seed)
{
	float h, s, c;
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	// Load constants and lig into shared memory.
	__syncthreads();
	sincosf(h, &s, &c);
}

kernel::kernel(const float* h_sf_e, const float* h_sf_d, const int h_sf_ns, const int h_sf_ne, const float* h_corner0, const float* h_corner1, const float* h_num_probes, const float h_granularity_inverse, const int num_mc_tasks, const int num_generations) : num_mc_tasks(num_mc_tasks), num_generations(num_generations)
{
	// Initialize scoring function.
	cudaMalloc(&d_sf_e, h_sf_ne);
	cudaMalloc(&d_sf_d, h_sf_ne);
	cudaMemcpy(d_sf_e, h_sf_e, sizeof(float) * h_sf_ne, cudaMemcpyHostToDevice);
	cudaMemcpy(d_sf_d, h_sf_d, sizeof(float) * h_sf_ne, cudaMemcpyHostToDevice);
	assert(sizeof(c_sf_e)  == sizeof(d_sf_e));
	assert(sizeof(c_sf_d)  == sizeof(d_sf_d));
	assert(sizeof(c_sf_ns) == sizeof(d_sf_ns));
	cudaMemcpyToSymbol(c_sf_e,  &d_sf_e,  sizeof(c_sf_e));
	cudaMemcpyToSymbol(c_sf_d,  &d_sf_d,  sizeof(c_sf_d));
	cudaMemcpyToSymbol(c_sf_ns, &h_sf_ns, sizeof(c_sf_ns));

	// Initialize receptor.
	assert(sizeof(c_corner0) == sizeof(float) * 3);
	assert(sizeof(c_corner1) == sizeof(float) * 3);
	assert(sizeof(c_num_probes) == sizeof(float) * 3);
	assert(sizeof(c_granularity_inverse) == sizeof(h_granularity_inverse));
	cudaMemcpyToSymbol(c_corner0, h_corner0, sizeof(c_corner0));
	cudaMemcpyToSymbol(c_corner1, h_corner1, sizeof(c_corner1));
	cudaMemcpyToSymbol(c_num_probes, h_num_probes, sizeof(c_num_probes));
	cudaMemcpyToSymbol(c_granularity_inverse, &h_granularity_inverse, sizeof(c_granularity_inverse));
	assert(sizeof(d_maps) == sizeof(float*) * sf_n);
	memset(d_maps, 0, sizeof(d_maps));
}

void kernel::update(const vector<vector<float> > h_maps, const size_t map_bytes, const vector<size_t>& xs)
{
	for (int i = 0; i < xs.size(); ++i)
	{
		const size_t t = xs[i];
		float* d_m;
		cudaMalloc(&d_m, map_bytes);
		cudaMemcpy(d_m, &h_maps[t].front(), map_bytes, cudaMemcpyHostToDevice);
		d_maps[t] = d_m;
	}
	assert(sizeof(c_maps) == sizeof(d_maps));
	cudaMemcpyToSymbol(c_maps, d_maps, sizeof(c_maps));
}

void kernel::launch(vector<float>& h_ex, const int* h_lig, const int nv, const int nf, const int na, const int np, const size_t* seed)
{
	// Copy ligand content from host memory to device memory.
	const size_t lig_bytes = sizeof(int) * (11 * nf + nf - 1 + 4 * na + 3 * np);
	float* d_lig;
	cudaMalloc(&d_lig, lig_bytes);
	cudaMemcpy(d_lig, &h_lig, lig_bytes, cudaMemcpyHostToDevice);

	// Allocate device memory for solutions.
	const size_t sln_bytes = sizeof(float) * (1 + (nv + 1) + nv + 3 * nf + 4 * nf + 3 * na + 3 * na + 3 * nf + 3 * nf) * num_mc_tasks;
	float* d_s0, d_s1, d_s2;
	cudaMalloc(&d_s0, sln_bytes);
	cudaMalloc(&d_s1, sln_bytes);
	cudaMalloc(&d_s2, sln_bytes);

	// Invoke CUDA kernel.
	bfgs<<<num_mc_tasks / 128, 128, lig_bytes>>>(d_s0, d_s1, d_s2, d_lig, nv, nf, na, np);

	// Copy e and x from device memory to host memory.
	const size_t ex_size = (1 + nv + 1) * num_mc_tasks;
	const size_t ex_bytes = sizeof(float) * ex_size;
	h_ex.resize(ex_size);
	cudaMemcpy(h_ex, d_s0, ex_bytes, cudaMemcpyDeviceToHost);

	// Free device memory.
	cudaFree(d_s0);
	cudaFree(d_s1);
	cudaFree(d_s2);
	cudaFree(d_lig);
}

kernel::~kernel()
{
	for (size_t t = 0; t < sf_n; ++t)
	{
		const float* d_m = d_maps[t];
		if (d_m) cudaFree(d_m);
	}
	cudaFree(d_sf_d);
	cudaFree(d_sf_e);
}

// -arch=sm_13 -maxrregcount=N -Xptxas=-v -ftz=true -prec-div=false -prec-sqrt=false -use_fast_math