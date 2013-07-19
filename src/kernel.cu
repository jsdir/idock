#include <assert.h>
#include <curand_kernel.h>
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
__constant__ int3 c_num_probes;
__constant__ float c_granularity_inverse;
__constant__ float* c_maps[sf_n];
__constant__ int c_num_generations;

extern __shared__ float shared[];

__device__  __noinline__// __forceinline__
bool evaluate(float* e, float* g, float* a, float* q, float* c, float* d, float* f, float* t, const float* x, const float eub)
{
	return true;
}

__global__
//__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
void bfgs(float* __restrict__ s0e, const int* lig, const int nv, const int nf, const int na, const int np, const unsigned long long seed)
{
	const int gid = blockIdx.x * blockDim.x + threadIdx.x;
	const int gds = blockDim.x * gridDim.x;
	const int nls = 5; // Number of line search trials for determining step size in BFGS
	const float eub = 40.0f * na; // A conformation will be droped if its free energy is not better than e_upper_bound.
	float* s0x = s0e + gds;
	float* s0g = s0x + (nv + 1) * gds;
	float* s0a = s0g + nv * gds;
	float* s0q = s0a + 3 * nf * gds;
	float* s0c = s0q + 4 * nf * gds;
	float* s0d = s0c + 3 * na * gds;
	float* s0f = s0d + 3 * na * gds;
	float* s0t = s0f + 3 * nf * gds;
	float* s1e = s0t + 3 * nf * gds;
	float* s1x = s1e + gds;
	float* s1g = s1x + (nv + 1) * gds;
	float* s1a = s1g + nv * gds;
	float* s1q = s1a + 3 * nf * gds;
	float* s1c = s1q + 4 * nf * gds;
	float* s1d = s1c + 3 * na * gds;
	float* s1f = s1d + 3 * na * gds;
	float* s1t = s1f + 3 * nf * gds;
	float* s2e = s1t + 3 * nf * gds;
	float* s2x = s2e + gds;
	float* s2g = s2x + (nv + 1) * gds;
	float* s2a = s2g + nv * gds;
	float* s2q = s2a + 3 * nf * gds;
	float* s2c = s2q + 4 * nf * gds;
	float* s2d = s2c + 3 * na * gds;
	float* s2f = s2d + 3 * na * gds;
	float* s2t = s2f + 3 * nf * gds;
	float* bfh = s2t + 3 * nf * gds;
	float* bfp = bfh + (nv*(nv+1)>>1) * gds;
	float* bfy = bfp + nv * gds;
	float* bfm = bfy + nv * gds;
	float rd0, rd1, rd2, rd3, rst;
	float sum, pg1, pga, pgc, alp, pg2, pr0, pr1, pr2, nrm, ang, sng, pq0, pq1, pq2, pq3, s1xq0, s1xq1, s1xq2, s1xq3, s2xq0, s2xq1, s2xq2, s2xq3, bpi;
	float yhy, yps, ryp, pco, bpj, bmj, ppj;
	int g, i, j, o0, o1, o2;
	curandState crs;

	// Load ligand into external shared memory.
	g = 11 * nf + nf - 1 + 4 * na + 3 * np;
	o0 = threadIdx.x;
	for (i = 0, j = (g - 1) / blockDim.x; i < j; ++i)
	{
		shared[o0] = lig[o0];
		o0 += blockDim.x;
	}
	if (o0 < g)
	{
		shared[o0] = lig[o0];
	}
	__syncthreads();

	curand_init(seed, gid, 0, &crs);
	curand_uniform(&crs);
}

kernel::kernel(const float* h_sf_e, const float* h_sf_d, const int h_sf_ns, const int h_sf_ne, const float* h_corner0, const float* h_corner1, const int* h_num_probes, const float h_granularity_inverse, const int num_mc_tasks, const int ng) : num_mc_tasks(num_mc_tasks), ng(ng)
{
	// Initialize scoring function.
	cudaMalloc(&d_sf_e, h_sf_ne);
	cudaMalloc(&d_sf_d, h_sf_ne);
	cudaMemcpy(d_sf_e, h_sf_e, sizeof(float) * h_sf_ne, cudaMemcpyHostToDevice);
	cudaMemcpy(d_sf_d, h_sf_d, sizeof(float) * h_sf_ne, cudaMemcpyHostToDevice);
	assert(sizeof(c_sf_e)  == sizeof(d_sf_e));
	assert(sizeof(c_sf_d)  == sizeof(d_sf_d));
	assert(sizeof(c_sf_ns) == sizeof(h_sf_ns));
	cudaMemcpyToSymbol(c_sf_e,  &d_sf_e,  sizeof(c_sf_e));
	cudaMemcpyToSymbol(c_sf_d,  &d_sf_d,  sizeof(c_sf_d));
	cudaMemcpyToSymbol(c_sf_ns, &h_sf_ns, sizeof(c_sf_ns));

	// Initialize receptor.
	assert(sizeof(c_corner0) == sizeof(float) * 3);
	assert(sizeof(c_corner1) == sizeof(float) * 3);
	assert(sizeof(c_num_probes) == sizeof(int) * 3);
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

void kernel::launch(vector<float>& h_ex, const vector<int>& h_lig, const int nv, const int nf, const int na, const int np, const unsigned long long seed)
{
	// Copy ligand content from host memory to device memory.
	const size_t lig_bytes = sizeof(int) * h_lig.size();
	int* d_lig;
	cudaMalloc(&d_lig, lig_bytes);
	cudaMemcpy(d_lig, &h_lig.front(), lig_bytes, cudaMemcpyHostToDevice);

	// Allocate device memory for variables. 3 * (nt + 1) is sufficient for t because the torques of inactive frames are always zero.
	const size_t var_bytes = sizeof(float) * (1 + (nv + 1) + nv + 3 * nf + 4 * nf + 3 * na + 3 * na + 3 * nf + 3 * nf) * num_mc_tasks * 3 + (nv * (nv + 1) >> 1) * num_mc_tasks + nv * num_mc_tasks * 3;
	float* d_s0;
	cudaMalloc(&d_s0, var_bytes);

	// Invoke CUDA kernel.
	bfgs<<<num_mc_tasks / 128, 128, lig_bytes>>>(d_s0, d_lig, nv, nf, na, np, seed);

	// Copy e and x from device memory to host memory.
	const size_t ex_size = (1 + nv + 1) * num_mc_tasks;
	const size_t ex_bytes = sizeof(float) * ex_size;
	h_ex.resize(ex_size);
	cudaMemcpy(&h_ex.front(), d_s0, ex_bytes, cudaMemcpyDeviceToHost);

	// Free device memory.
	cudaFree(d_s0);
	cudaFree(d_lig);
}

kernel::~kernel()
{
	for (size_t t = 0; t < sf_n; ++t)
	{
		float* const d_m = d_maps[t];
		if (d_m) cudaFree(d_m);
	}
	cudaFree(d_sf_d);
	cudaFree(d_sf_e);
}

// -arch=sm_13 -maxrregcount=N -Xptxas=-v -ftz=true -prec-div=false -prec-sqrt=false -use_fast_math