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
void bfgs(float* __restrict__ s0e, float* __restrict__ s1e, float* __restrict__ s2e, const int seed)
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

void kernel::update(const float* h_maps[], const size_t num_probes_product, const size_t* xs, const size_t n)
{
	for (int i = 0; i < n; ++i)
	{
		const size_t t = xs[i];
		float* d_m;
		cudaMalloc(&d_m, sizeof(float) * num_probes_product);
		cudaMemcpy(d_m, &h_maps[t], sizeof(float) * num_probes_product, cudaMemcpyHostToDevice);
		d_maps[t] = d_m;
	}
	cudaMemcpyToSymbol(c_maps, &d_maps.front(), sizeof(c_maps));
}

void kernel::launch(float* h_s0, const float* h_lig, const int n, const size_t oz, const size_t og, const size_t* seed)
{
//	cudaMemcpy(d_lig, &h_lig, sizeof(float) * n, cudaMemcpyHostToDevice);

	const size_t solution_size = oz * num_mc_tasks;
	float* d_s0, d_s1, d_s2;
	cudaMalloc(&d_s0, solution_size);
	cudaMalloc(&d_s1, solution_size);
	cudaMalloc(&d_s2, solution_size);

	bfgs<<<num_mc_tasks / 128, 128, sizeof(float) * n>>>(d_s0, d_s1, d_s2);

	// Transfer e and x only.
	const size_t result_size = og * num_mc_tasks;
	h_s0.resize(result_size);
	cudaMemcpy(h_s0, d_s0, result_size, cudaMemcpyDeviceToHost);

	// Free device memory.
	cudaFree(d_s0);
	cudaFree(d_s1);
	cudaFree(d_s2);
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