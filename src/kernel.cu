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
__constant__ int c_ng;

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

	// Randomize s0x.
	curand_init(seed, gid, 0, &crs);
	rd0 = curand_uniform(&crs) * 2 - 1;
	s0x[o0  = gid] = 0.5f * ((1 + rd0) * c_corner1[0] + (1 - rd0) * c_corner0[0]);
	rd0 = curand_uniform(&crs) * 2 - 1;
	s0x[o0 += gds] = 0.5f * ((1 + rd0) * c_corner1[1] + (1 - rd0) * c_corner0[1]);
	rd0 = curand_uniform(&crs) * 2 - 1;
	s0x[o0 += gds] = 0.5f * ((1 + rd0) * c_corner1[2] + (1 - rd0) * c_corner0[2]);
	rd0 = curand_uniform(&crs) * 2 - 1;
	rd1 = curand_uniform(&crs) * 2 - 1;
	rd2 = curand_uniform(&crs) * 2 - 1;
	rd3 = curand_uniform(&crs) * 2 - 1;
	rst = 1.0f / sqrt(rd0*rd0 + rd1*rd1 + rd2*rd2 + rd3*rd3);
//	rst = rsqrtf(rd0*rd0 + rd1*rd1 + rd2*rd2 + rd3*rd3);
	s0x[o0 += gds] = rd0 * rst;
	s0x[o0 += gds] = rd1 * rst;
	s0x[o0 += gds] = rd2 * rst;
	s0x[o0 += gds] = rd3 * rst;
	for (i = 6; i < nv; ++i)
	{
		s0x[o0 += gds] = curand_uniform(&crs) * 2 - 1;
	}
	evaluate(s0e, s0g, s0a, s0q, s0c, s0d, s0f, s0t, s0x, eub);

	// Repeat for a number of generations.
	for (g = 0; g < ng; ++g)
	{
		// Mutate s0x into s1x
		o0  = gid;
		s1x[o0] = s0x[o0] + curand_uniform(&crs) * 2 - 1;
		o0 += gds;
		s1x[o0] = s0x[o0] + curand_uniform(&crs) * 2 - 1;
		o0 += gds;
		s1x[o0] = s0x[o0] + curand_uniform(&crs) * 2 - 1;
//		for (i = 3; i < nv + 1; ++i)
		for (i = 2 - nv; i < 0; ++i)
		{
			o0 += gds;
			s1x[o0] = s0x[o0];
		}
		evaluate(s1e, s1g, s1a, s1q, s1c, s1d, s1f, s1t, s1x, eub);

		// Initialize the inverse Hessian matrix to identity matrix.
		// An easier option that works fine in practice is to use a scalar multiple of the identity matrix,
		// where the scaling factor is chosen to be in the range of the eigenvalues of the true Hessian.
		// See N&R for a recipe to find this initializer.
		bfh[o0 = gid] = 1.0f;
		for (j = 1; j < nv; ++j)
		{
			for (i = 0; i < j; ++i)
			{
				bfh[o0 += gds] = 0.0f;
			}
			bfh[o0 += gds] = 1.0f;
		}

		// Use BFGS to optimize the mutated conformation s1x into local optimum s2x.
		// http://en.wikipedia.org/wiki/BFGS_method
		// http://en.wikipedia.org/wiki/Quasi-Newton_method
		// The loop breaks when no appropriate alpha can be found.
		while (true)
		{
			// Calculate p = -h * g, where p is for descent direction, h for Hessian, and g for gradient.
			sum = bfh[o1 = gid] * s1g[o0 = gid];
			for (i = 1; i < nv; ++i)
			{
				sum += bfh[o1 += i * gds] * s1g[o0 += gds];
			}
			bfp[o2 = gid] = -sum;
			for (j = 1; j < nv; ++j)
			{
				sum = bfh[o1 = (j*(j+1)>>1) * gds + gid] * s1g[o0 = gid];
				for (i = 1; i < nv; ++i)
				{
					sum += bfh[o1 += i > j ? i * gds : gds] * s1g[o0 += gds];
				}
				bfp[o2 += gds] = -sum;
			}

			// Calculate pg = p * g = -h * g^2 < 0
			o0 = gid;
			pg1 = bfp[o0] * s1g[o0];
			for (i = 1; i < nv; ++i)
			{
				o0 += gds;
				pg1 += bfp[o0] * s1g[o0];
			}
			pga = 0.0001f * pg1;
			pgc = 0.9f * pg1;

			// Perform a line search to find an appropriate alpha.
			// Try different alpha values for nls times.
			// alpha starts with 1, and shrinks to 0.1 of itself iteration by iteration.
			alp = 1.0f;
			for (j = 0; j < nls; ++j)
			{
				// Calculate x2 = x1 + a * p.
				o0  = gid;
				s2x[o0] = s1x[o0] + alp * bfp[o0];
				o0 += gds;
				s2x[o0] = s1x[o0] + alp * bfp[o0];
				o0 += gds;
				s2x[o0] = s1x[o0] + alp * bfp[o0];
				o0 += gds;
				s1xq0 = s1x[o0];
				pr0 = bfp[o0];
				o0 += gds;
				s1xq1 = s1x[o0];
				pr1 = bfp[o0];
				o0 += gds;
				s1xq2 = s1x[o0];
				pr2 = bfp[o0];
				o0 += gds;
				s1xq3 = s1x[o0];
				assert(fabs(s1xq0*s1xq0 + s1xq1*s1xq1 + s1xq2*s1xq2 + s1xq3*s1xq3 - 1.0f) < 1e-3f);
				nrm = sqrt(pr0*pr0 + pr1*pr1 + pr2*pr2);
				ang = 0.5f * alp * nrm;
				sng = sinf(ang) / nrm;
				pq0 = cosf(ang);
//				sincosf(ang, &sng, &pq0);
//				sincospif(ang, &sng, &pq0);
//				sng /= nrm;
				pq1 = sng * pr0;
				pq2 = sng * pr1;
				pq3 = sng * pr2;
				assert(fabs(pq0*pq0 + pq1*pq1 + pq2*pq2 + pq3*pq3 - 1.0f) < 1e-3f);
				s2xq0 = pq0 * s1xq0 - pq1 * s1xq1 - pq2 * s1xq2 - pq3 * s1xq3;
				s2xq1 = pq0 * s1xq1 + pq1 * s1xq0 + pq2 * s1xq3 - pq3 * s1xq2;
				s2xq2 = pq0 * s1xq2 - pq1 * s1xq3 + pq2 * s1xq0 + pq3 * s1xq1;
				s2xq3 = pq0 * s1xq3 + pq1 * s1xq2 - pq2 * s1xq1 + pq3 * s1xq0;
				assert(fabs(s2xq0*s2xq0 + s2xq1*s2xq1 + s2xq2*s2xq2 + s2xq3*s2xq3 - 1.0f) < 1e-3f);
				s2x[o0 -= 3 * gds] = s2xq0;
				s2x[o0 += gds] = s2xq1;
				s2x[o0 += gds] = s2xq2;
				s2x[o0 += gds] = s2xq3;
				for (i = 6; i < nv; ++i)
				{
					bpi = bfp[o0];
					o0 += gds;
					s2x[o0] = s1x[o0] + alp * bpi;
				}

				// Evaluate x2, subject to Wolfe conditions http://en.wikipedia.org/wiki/Wolfe_conditions
				// 1) Armijo rule ensures that the step length alpha decreases f sufficiently.
				// 2) The curvature condition ensures that the slope has been reduced sufficiently.
				if (evaluate(s2e, s2g, s2a, s2q, s2c, s2d, s2f, s2t, s2x, s1e[gid] + alp * pga))
				{
					o0 = gid;
					pg2 = bfp[o0] * s2g[o0];
					for (i = 1; i < nv; ++i)
					{
						o0 += gds;
						pg2 += bfp[o0] * s2g[o0];
					}
					if (pg2 >= pgc) break;
				}

				alp *= 0.1f;
			}

			// If no appropriate alpha can be found, exit the BFGS loop.
			if (j == nls) break;

			// Calculate y = g2 - g1.
			o0 = gid;
			bfy[o0] = s2g[o0] - s1g[o0];
			for (i = 1; i < nv; ++i)
			{
				o0 += gds;
				bfy[o0] = s2g[o0] - s1g[o0];
			}

			// Calculate m = -h * y.
			sum = bfh[o1 = gid] * bfy[o0 = gid];
			for (i = 1; i < nv; ++i)
			{
				sum += bfh[o1 += i * gds] * bfy[o0 += gds];
			}
			bfm[o2 = gid] = -sum;
			for (j = 1; j < nv; ++j)
			{
				sum = bfh[o1 = (j*(j+1)>>1) * gds + gid] * bfy[o0 = gid];
				for (i = 1; i < nv; ++i)
				{
					sum += bfh[o1 += i > j ? i * gds : gds] * bfy[o0 += gds];
				}
				bfm[o2 += gds] = -sum;
			}

			// Calculate yhy = -y * m = -y * (-h * y) = y * h * y.
			o0 = gid;
			yhy = -bfy[o0] * bfm[o0];
			for (i = 1; i < nv; ++i)
			{
				o0 += gds;
				yhy -= bfy[o0] * bfm[o0];
			}

			// Calculate yps = y * p.
			o0 = gid;
			yps = bfy[o0] * bfp[o0];
			for (i = 1; i < nv; ++i)
			{
				o0 += gds;
				yps += bfy[o0] * bfp[o0];
			}

			// Update Hessian matrix h.
			ryp = 1.0f / yps;
			pco = ryp * (ryp * yhy + alp);
			o2 = gid;
			for (j = 0; j < nv; ++j)
			{
				bpj = bfp[o2];
				bmj = bfm[o2];
				ppj = pco * bpj;
				bfh[o1 = (j*(j+3)>>1) * gds + gid] += (ryp * 2 * bmj + ppj) * bpj;
				for (i = j + 1; i < nv; ++i)
				{
					o0 = i * gds + gid;
					bpi = bfp[o0];
					bfh[o1 += i * gds] += ryp * (bmj * bpi + bfm[o0] * bpj) + ppj * bpi;
				}
				o2 += gds;
			}

			// Move to the next iteration, i.e. e1 = e2, x1 = x2, g1 = g2.
			o0 = gid;
			s1e[o0] = s2e[o0];
//			for (i = 1; i < 2 * (nv + 1); ++i)
			for (i = -1 - 2 * nv; i < 0; ++i)
			{
				o0 += gds;
				s1e[o0] = s2e[o0];
			}
		}

		// Accept x1 according to Metropolis criteria.
		if (s1e[gid] < s0e[gid])
		{
			o0 = gid;
			s0e[o0] = s1e[o0];
//			for (i = 1; i < nv + 2; ++i)
			for (i = -1 - nv; i < 0; ++i)
			{
				o0 += gds;
				s0e[o0] = s1e[o0];
			}
		}
	}
}

kernel::kernel(const float* h_sf_e, const float* h_sf_d, const int h_sf_ns, const int h_sf_ne, const float* h_corner0, const float* h_corner1, const int* h_num_probes, const float h_granularity_inverse, const int num_mc_tasks, const int h_ng) : num_mc_tasks(num_mc_tasks)
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
	assert(sizeof(c_ng) == sizeof(h_ng));
	cudaMemcpyToSymbol(c_corner0, h_corner0, sizeof(c_corner0));
	cudaMemcpyToSymbol(c_corner1, h_corner1, sizeof(c_corner1));
	cudaMemcpyToSymbol(c_num_probes, h_num_probes, sizeof(c_num_probes));
	cudaMemcpyToSymbol(c_granularity_inverse, &h_granularity_inverse, sizeof(c_granularity_inverse));
	cudaMemcpyToSymbol(c_ng, &h_ng, sizeof(c_ng));
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