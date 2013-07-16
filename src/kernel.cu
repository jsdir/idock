#include <assert.h>
#include "kernel.hpp"

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#undef  assert
#define assert(arg)
#endif

__constant__ float* c_sfe;
__constant__ float* c_sfd;
__constant__ float3 c_corner0;
__constant__ float3 c_corner1;
__constant__ float3 c_num_probes;
__constant__ float c_granularity_inverse;
__constant__ float* c_maps[scoring_function::n];
__constant__ int c_num_generations;

extern __shared__ float shared[];

__device__  __noinline__// __forceinline__
bool evaluate(const float* x, float* e, float* g, float* a, float* q, float* c, float* d, float* f, float* t, const float e_upper_bound)
{
	return true;
}

__global__
//__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
void bfgs(float* __restrict__ s0e, float* __restrict__ s1e, float* __restrict__ s2e, const size_t seed)
{
	float h, s, c;
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	// Load constants into shared memory.
	__syncthreads();
	sincosf(h, &s, &c);
}

kernel::kernel(const scoring_function& sf, const receptor& rec, const size_t num_mc_tasks, const size_t num_generations) : num_mc_tasks(num_mc_tasks), num_generations(num_generations)
{
	cudaMalloc(&d_sfe, sf.ne);
	cudaMalloc(&d_sfd, sf.ne);
	cudaMemcpyToSymbol(c_sfe, &d_sfe, sizeof(c_sfe));
	cudaMemcpyToSymbol(c_sfd, &d_sfd, sizeof(c_sfd));
	cudaMemcpyToSymbol(c_corner0, &rec.corner0.front(), sizeof(c_corner0));
	cudaMemcpyToSymbol(c_corner1, &rec.corner1.front(), sizeof(c_corner1));
	cudaMemcpyToSymbol(c_num_probes, &rec.num_probes.front(), sizeof(c_num_probes));
	cudaMemcpyToSymbol(c_granularity_inverse, &rec.granularity_inverse, sizeof(c_granularity_inverse));
	memset(d_maps, 0, sizeof(d_maps));
}

void kernel::update(const receptor& rec, const vector<size_t>& xs)
{
	for (size_t i = 0; i < xs.size(); ++i)
	{
		const size_t t = xs[i];
		float* d_m;
		cudaMalloc(&d_m, sizeof(float) * rec.num_probes_product);
		cudaMemcpy(d_m, &rec.maps[t].front(), sizeof(float) * rec.num_probes_product, cudaMemcpyHostToDevice);
		d_maps[t] = d_m;
	}
	cudaMemcpyToSymbol(c_maps, &d_maps.front(), sizeof(c_maps));
}

void kernel::launch(vector<float>& h_s0, const ligand& lig)
{
	const size_t solution_size = lig.oz * num_mc_tasks;
	float* d_s0, d_s1, d_s2;
	cudaMalloc(&d_s0, solution_size);
	cudaMalloc(&d_s1, solution_size);
	cudaMalloc(&d_s2, solution_size);
	bfgs<<<num_mc_tasks / 128, 128/*, sm_size*/>>>(d_s0, d_s1, d_s2);
	h_s0.resize(solution_size); // Only copy e and x.
	cudaMemcpy(&h_s0.front(), d_s0, solution_size, cudaMemcpyDeviceToHost);
	cudaFree(d_s0);
	cudaFree(d_s1);
	cudaFree(d_s2);
}

kernel::~kernel()
{
	for (size_t t = 0; t < scoring_function::n; ++t)
	{
		const float* d_m = d_maps[t];
		if (d_m) cudaFree(d_m);
	}
	cudaFree(d_sfd);
	cudaFree(d_sfe);
}

// -arch=sm_13 -maxrregcount=N -Xptxas=-v -ftz=true -prec-div=false -prec-sqrt=false -use_fast_math