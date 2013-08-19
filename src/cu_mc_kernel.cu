#include <cassert>
#include <curand_kernel.h>
#include <helper_cuda.h>
#include "cu_mc_kernel.hpp"

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#undef  assert
#define assert(arg)
#endif

__constant__ const float* c_sf_e;
__constant__ const float* c_sf_d;
__constant__ int c_sf_ns;
__constant__ float3 c_corner0;
__constant__ float3 c_corner1;
__constant__ int3 c_num_probes;
__constant__ float c_granularity_inverse;
__constant__ const float* c_maps[sf_n];
__constant__ int c_ng;
__constant__ unsigned long c_seed;

extern __shared__ int shared[];

// TODO: Use Restrict Qualifier for Kernel Arguments
__device__  __noinline__// __forceinline__
bool evaluate(float* e, float* g, float* a, float* q, float* c, float* d, float* f, float* t, const float* x, const int nf, const int na, const int np, const float eub)
{
	const int gid = blockIdx.x * blockDim.x + threadIdx.x;
	const int gds = blockDim.x * gridDim.x;
	const int gd3 = 3 * gds;
	const int gd4 = 4 * gds;
	const int* const act = shared;
	const int* const beg = &act[nf];
	const int* const end = &beg[nf];
	const int* const nbr = &end[nf];
	const int* const prn = &nbr[nf];
	const float* const yy0 = (float*)&prn[nf];
	const float* const yy1 = &yy0[nf];
	const float* const yy2 = &yy1[nf];
	const float* const xy0 = &yy2[nf];
	const float* const xy1 = &xy0[nf];
	const float* const xy2 = &xy1[nf];
	const int* const brs = (int*)&xy2[nf];
	const float* const cr0 = (float*)&brs[nf - 1];
	const float* const cr1 = &cr0[na];
	const float* const cr2 = &cr1[na];
	const int* const xst = (int*)&cr2[na];
	const int* const ip0 = &xst[na];
	const int* const ip1 = &ip0[np];
	const int* const ipp = &ip1[np];

	float y, y0, y1, y2, v0, v1, v2, c0, c1, c2, e000, e100, e010, e001, a0, a1, a2, ang, sng, r0, r1, r2, r3, vs, dr, f0, f1, f2, t0, t1, t2, d0, d1, d2;
	float q0, q1, q2, q3, q00, q01, q02, q03, q11, q12, q13, q22, q23, q33, m0, m1, m2, m3, m4, m5, m6, m7, m8;
	int i, j, k, b, w, i0, i1, i2, k0, k1, k2, z;
	const float* map;

	// Apply position, orientation and torsions.
	c[i  = gid] = x[k  = gid];
	c[i += gds] = x[k += gds];
	c[i += gds] = x[k += gds];
	q[i  = gid] = x[k += gds];
	q[i += gds] = x[k += gds];
	q[i += gds] = x[k += gds];
	q[i += gds] = x[k += gds];
	y = 0.0f;
	for (k = 0, b = 0, w = 6 * gds + gid; k < nf; ++k)
	{
		// Load rotorY from memory into registers.
		y0 = c[i0  = beg[k] * gd3 + gid];
		y1 = c[i0 += gds];
		y2 = c[i0 += gds];

		// Translate orientation of active frames from quaternion into 3x3 matrix.
		if (act[k])
		{
			q0 = q[k0  = k * gd4 + gid];
			q1 = q[k0 += gds];
			q2 = q[k0 += gds];
			q3 = q[k0 += gds];
			assert(fabs(q0*q0 + q1*q1 + q2*q2 + q3*q3 - 1.0f) < 1e-3f);
			q00 = q0 * q0;
			q01 = q0 * q1;
			q02 = q0 * q2;
			q03 = q0 * q3;
			q11 = q1 * q1;
			q12 = q1 * q2;
			q13 = q1 * q3;
			q22 = q2 * q2;
			q23 = q2 * q3;
			q33 = q3 * q3;
			m0 = q00 + q11 - q22 - q33;
			m1 = 2 * (q12 - q03);
			m2 = 2 * (q02 + q13);
			m3 = 2 * (q03 + q12);
			m4 = q00 - q11 + q22 - q33;
			m5 = 2 * (q23 - q01);
			m6 = 2 * (q13 - q02);
			m7 = 2 * (q01 + q23);
			m8 = q00 - q11 - q22 + q33;
		}

		// Evaluate c and d of frame atoms. Aggregate e into y.
		for (i = beg[k], z = end[k]; i < z; ++i)
		{
			i0 = i * gd3 + gid;
			i1 = i0 + gds;
			i2 = i1 + gds;

			// The first atom of a frame is assumed to be its rotor Y.
			// TODO: avoid use of branching in short computations
			if (i == beg[k])
			{
				c0 = y0;
				c1 = y1;
				c2 = y2;
			}
			else
			{
				// Calculate coordinate from transformation matrix and offset.
				v0 = cr0[i];
				v1 = cr1[i];
				v2 = cr2[i];
				c0 = y0 + m0 * v0 + m1 * v1 + m2 * v2;
				c1 = y1 + m3 * v0 + m4 * v1 + m5 * v2;
				c2 = y2 + m6 * v0 + m7 * v1 + m8 * v2;

				// Store coordinate from registers into memory.
				c[i0] = c0;
				c[i1] = c1;
				c[i2] = c2;
			}

			// Penalize out-of-box case.
			if (c0 < c_corner0.x || c_corner1.x <= c0 || c1 < c_corner0.y || c_corner1.y <= c1 || c2 < c_corner0.z || c_corner1.z <= c2)
			{
				y += 10.0f;
				d[i0] = 0.0f;
				d[i1] = 0.0f;
				d[i2] = 0.0f;
				continue;
			}

			// Find the index of the current coordinate
			k0 = static_cast<int>((c0 - c_corner0.x) * c_granularity_inverse);
			k1 = static_cast<int>((c1 - c_corner0.y) * c_granularity_inverse);
			k2 = static_cast<int>((c2 - c_corner0.z) * c_granularity_inverse);
			assert(k0 + 1 < c_num_probes.x);
			assert(k1 + 1 < c_num_probes.y);
			assert(k2 + 1 < c_num_probes.z);
			k0 = c_num_probes.x * (c_num_probes.y * k2 + k1) + k0;

			// Retrieve the grid map and lookup the value
			map = c_maps[xst[i]];
			e000 = map[k0];
			e100 = map[k0 + 1];
			e010 = map[k0 + c_num_probes.x];
			e001 = map[k0 + c_num_probes.x * c_num_probes.y];
			y += e000;
			d[i0] = (e100 - e000) * c_granularity_inverse;
			d[i1] = (e010 - e000) * c_granularity_inverse;
			d[i2] = (e001 - e000) * c_granularity_inverse;
		}
		for (j = 0, z = nbr[k]; j < z; ++j)
		{
			i = brs[b++];
			i0 = beg[i] * gd3 + gid;
			i1 = i0 + gds;
			i2 = i1 + gds;
			c[i0] = y0 + m0 * yy0[i] + m1 * yy1[i] + m2 * yy2[i];
			c[i1] = y1 + m3 * yy0[i] + m4 * yy1[i] + m5 * yy2[i];
			c[i2] = y2 + m6 * yy0[i] + m7 * yy1[i] + m8 * yy2[i];

			// Skip inactive BRANCH frame
			if (!act[i]) continue;

			// Update a of BRANCH frame
			a0 = m0 * xy0[i] + m1 * xy1[i] + m2 * xy2[i];
			a1 = m3 * xy0[i] + m4 * xy1[i] + m5 * xy2[i];
			a2 = m6 * xy0[i] + m7 * xy1[i] + m8 * xy2[i];
			assert(fabs(a0*a0 + a1*a1 + a2*a2 - 1.0f) < 1e-3f);
			a[k0  = i * gd3 + gid] = a0;
			a[k0 += gds] = a1;
			a[k0 += gds] = a2;

			// Update q of BRANCH frame
			ang = x[w += gds] * 0.5f;
//			sng = sinf(ang);
//			r0 = cosf(ang);
			sincosf(ang, &sng, &r0);
//			sincospif(ang, &sng, &r0);
			r1 = sng * a0;
			r2 = sng * a1;
			r3 = sng * a2;
			q00 = r0 * q0 - r1 * q1 - r2 * q2 - r3 * q3;
			q01 = r0 * q1 + r1 * q0 + r2 * q3 - r3 * q2;
			q02 = r0 * q2 - r1 * q3 + r2 * q0 + r3 * q1;
			q03 = r0 * q3 + r1 * q2 - r2 * q1 + r3 * q0;
			assert(fabs(q00*q00 + q01*q01 + q02*q02 + q03*q03 - 1.0f) < 1e-3f);
			q[k0  = i * gd4 + gid] = q00;
			q[k0 += gds] = q01;
			q[k0 += gds] = q02;
			q[k0 += gds] = q03;
		}
	}
	assert(b == nf - 1);
//	assert(w == nv * gds + gid);
	assert(k == nf);

	// Calculate intra-ligand free energy.
	for (i = 0; i < np; ++i)
	{
		i0 = ip0[i] * gd3 + gid;
		i1 = i0 + gds;
		i2 = i1 + gds;
		k0 = ip1[i] * gd3 + gid;
		k1 = k0 + gds;
		k2 = k1 + gds;
		v0 = c[k0] - c[i0];
		v1 = c[k1] - c[i1];
		v2 = c[k2] - c[i2];
		vs = v0*v0 + v1*v1 + v2*v2;
		if (vs < 64.0)
		{
			j = ipp[i] + static_cast<int>(c_sf_ns * vs);
			y += c_sf_e[j];
			dr = c_sf_d[j];
			d0 = dr * v0;
			d1 = dr * v1;
			d2 = dr * v2;
			d[i0] -= d0;
			d[i1] -= d1;
			d[i2] -= d2;
			d[k0] += d0;
			d[k1] += d1;
			d[k2] += d2;
		}
	}

	// If the free energy is no better than the upper bound, refuse this conformation.
	if (y >= eub) return false;

	// Store e from register into memory.
	e[gid] = y;

	// Calculate and aggregate the force and torque of BRANCH frames to their parent frame.
	f[k0 = gid] = 0.0f;
	t[k0] = 0.0f;
	for (i = 1, z = 3 * nf; i < z; ++i)
	{
		f[k0 += gds] = 0.0f;
		t[k0] = 0.0f;
	}
//	assert(w == nv * gds + gid);
	assert(k == nf);
	while (k)
	{
		--k;

		// Load f, t and rotorY from memory into register
		k0 = k * gd3 + gid;
		k1 = k0 + gds;
		k2 = k1 + gds;
		f0 = f[k0];
		f1 = f[k1];
		f2 = f[k2];
		t0 = t[k0];
		t1 = t[k1];
		t2 = t[k2];
		y0 = c[i0  = beg[k] * gd3 + gid];
		y1 = c[i0 += gds];
		y2 = c[i0 += gds];

		// Aggregate frame atoms.
		for (i = beg[k], z = end[k]; i < z; ++i)
		{
			i0 = i * gd3 + gid;
			i1 = i0 + gds;
			i2 = i1 + gds;
			d0 = d[i0];
			d1 = d[i1];
			d2 = d[i2];

			// The derivatives with respect to the position, orientation, and torsions
			// would be the negative total force acting on the ligand,
			// the negative total torque, and the negative torque projections, respectively,
			// where the projections refer to the torque applied to the branch moved by the torsion,
			// projected on its rotation axi
			f0 += d0;
			f1 += d1;
			f2 += d2;
			if (i == beg[k]) continue;

			v0 = c[i0] - y0;
			v1 = c[i1] - y1;
			v2 = c[i2] - y2;
			t0 += v1 * d2 - v2 * d1;
			t1 += v2 * d0 - v0 * d2;
			t2 += v0 * d1 - v1 * d0;
		}

		if (k)
		{
			// Save the aggregated torque of active BRANCH frames to g.
			if (act[k])
			{
				g[w -= gds] = t0 * a[k0] + t1 * a[k1] + t2 * a[k2]; // dot product
			}

			// Aggregate the force and torque of current frame to its parent frame.
			k0 = prn[k] * gd3 + gid;
			k1 = k0 + gds;
			k2 = k1 + gds;
			f[k0] += f0;
			f[k1] += f1;
			f[k2] += f2;
			v0 = y0 - c[i0  = beg[prn[k]] * gd3 + gid];
			v1 = y1 - c[i0 += gds];
			v2 = y2 - c[i0 += gds];
			t[k0] += t0 + v1 * f2 - v2 * f1;
			t[k1] += t1 + v2 * f0 - v0 * f2;
			t[k2] += t2 + v0 * f1 - v1 * f0;
		}
	}
	assert(w == 6 * gds + gid);

	// Save the aggregated force and torque of ROOT frame to g.
	g[i0  = gid] = f0;
	g[i0 += gds] = f1;
	g[i0 += gds] = f2;
	g[i0 += gds] = t0;
	g[i0 += gds] = t1;
	g[i0 += gds] = t2;
	return true;
}

__global__
//__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor) // .maxntid nx .minnctapersm ncta
void mc(float* __restrict__ s0e, const int* __restrict__ lig, const int nv, const int nf, const int na, const int np)
{
	const int gid = blockIdx.x * blockDim.x + threadIdx.x;
	const int gds = blockDim.x * gridDim.x;
	const int nls = 5; // Number of line search trials for determining step size in BFGS
	const float eub = 40.0f * na; // A conformation will be droped if its free energy is not better than e_upper_bound.
	float* const s0x = &s0e[gds];
	float* const s0g = &s0x[(nv + 1) * gds];
	float* const s0a = &s0g[nv * gds];
	float* const s0q = &s0a[3 * nf * gds];
	float* const s0c = &s0q[4 * nf * gds];
	float* const s0d = &s0c[3 * na * gds];
	float* const s0f = &s0d[3 * na * gds];
	float* const s0t = &s0f[3 * nf * gds];
	float* const s1e = &s0t[3 * nf * gds];
	float* const s1x = &s1e[gds];
	float* const s1g = &s1x[(nv + 1) * gds];
	float* const s1a = &s1g[nv * gds];
	float* const s1q = &s1a[3 * nf * gds];
	float* const s1c = &s1q[4 * nf * gds];
	float* const s1d = &s1c[3 * na * gds];
	float* const s1f = &s1d[3 * na * gds];
	float* const s1t = &s1f[3 * nf * gds];
	float* const s2e = &s1t[3 * nf * gds];
	float* const s2x = &s2e[gds];
	float* const s2g = &s2x[(nv + 1) * gds];
	float* const s2a = &s2g[nv * gds];
	float* const s2q = &s2a[3 * nf * gds];
	float* const s2c = &s2q[4 * nf * gds];
	float* const s2d = &s2c[3 * na * gds];
	float* const s2f = &s2d[3 * na * gds];
	float* const s2t = &s2f[3 * nf * gds];
	float* const bfh = &s2t[3 * nf * gds];
	float* const bfp = &bfh[(nv*(nv+1)>>1) * gds];
	float* const bfy = &bfp[nv * gds];
	float* const bfm = &bfy[nv * gds];
	float rd0, rd1, rd2, rd3, rst;
	float sum, pg1, pga, pgc, alp, pg2, pr0, pr1, pr2, nrm, ang, sng, pq0, pq1, pq2, pq3, s1xq0, s1xq1, s1xq2, s1xq3, s2xq0, s2xq1, s2xq2, s2xq3, bpi;
	float yhy, yps, ryp, pco, bpj, bmj, ppj;
	int g, i, j, o0, o1, o2;
	curandState crs;

	// Load ligand into external shared memory.
	// TODO: try not caching in shared memory.
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
	curand_init(c_seed, gid, 0, &crs);
	rd0 = curand_uniform(&crs);
	s0x[o0  = gid] = rd0 * c_corner1.x + (1 - rd0) * c_corner0.x;
	rd0 = curand_uniform(&crs);
	s0x[o0 += gds] = rd0 * c_corner1.y + (1 - rd0) * c_corner0.y;
	rd0 = curand_uniform(&crs);
	s0x[o0 += gds] = rd0 * c_corner1.z + (1 - rd0) * c_corner0.z;
	rd0 = curand_uniform(&crs);
	rd1 = curand_uniform(&crs);
	rd2 = curand_uniform(&crs);
	rd3 = curand_uniform(&crs);
	rst = rsqrtf(rd0*rd0 + rd1*rd1 + rd2*rd2 + rd3*rd3);
	s0x[o0 += gds] = rd0 * rst;
	s0x[o0 += gds] = rd1 * rst;
	s0x[o0 += gds] = rd2 * rst;
	s0x[o0 += gds] = rd3 * rst;
	for (i = 6; i < nv; ++i)
	{
		s0x[o0 += gds] = curand_uniform(&crs);
	}
/*
	s0x[o0  = gid] =  49.799f;
	s0x[o0 += gds] = -31.025f;
	s0x[o0 += gds] =  35.312f;
	s0x[o0 += gds] = 1.0f;
	s0x[o0 += gds] = 0.0f;
	s0x[o0 += gds] = 0.0f;
	s0x[o0 += gds] = 0.0f;
	for (i = 6; i < nv; ++i)
	{
		s0x[o0 += gds] = 0.0f;
	}
*/
	evaluate(s0e, s0g, s0a, s0q, s0c, s0d, s0f, s0t, s0x, nf, na, np, eub);

	// Mutate s0x into s1x
	o0  = gid;
	s1x[o0] = s0x[o0] + curand_uniform(&crs);
	o0 += gds;
	s1x[o0] = s0x[o0] + curand_uniform(&crs);
	o0 += gds;
	s1x[o0] = s0x[o0] + curand_uniform(&crs);
//	for (i = 3; i < nv + 1; ++i)
	for (i = 2 - nv; i < 0; ++i)
	{
		o0 += gds;
		s1x[o0] = s0x[o0];
	}
	evaluate(s1e, s1g, s1a, s1q, s1c, s1d, s1f, s1t, s1x, nf, na, np, eub);

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

	// Repeat for a number of generations.
	for (g = 0; g < c_ng; ++g)
	{
		// Use BFGS to optimize the mutated conformation s1x into local optimum s2x.
		// http://en.wikipedia.org/wiki/BFGS_method
		// http://en.wikipedia.org/wiki/Quasi-Newton_method

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
//			sng = sinf(ang) / nrm;
//			pq0 = cosf(ang);
			sincosf(ang, &sng, &pq0);
//			sincospif(ang, &sng, &pq0);
			sng /= nrm;
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
			if (evaluate(s2e, s2g, s2a, s2q, s2c, s2d, s2f, s2t, s2x, nf, na, np, s1e[gid] + alp * pga))
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

		// If no appropriate alpha can be found, restart the BFGS loop.
		if (j == nls)
		{
			// Accept x1 according to Metropolis criteria.
			if (s1e[gid] < s0e[gid])
			{
				o0 = gid;
				s0e[o0] = s1e[o0];
//				for (i = 1; i < nv + 2; ++i)
				for (i = -1 - nv; i < 0; ++i)
				{
					o0 += gds;
					s0e[o0] = s1e[o0];
				}
			}

			// Mutate s0x into s1x
			o0  = gid;
			s1x[o0] = s0x[o0] + curand_uniform(&crs);
			o0 += gds;
			s1x[o0] = s0x[o0] + curand_uniform(&crs);
			o0 += gds;
			s1x[o0] = s0x[o0] + curand_uniform(&crs);
//			for (i = 3; i < nv + 1; ++i)
			for (i = 2 - nv; i < 0; ++i)
			{
				o0 += gds;
				s1x[o0] = s0x[o0];
			}
			evaluate(s1e, s1g, s1a, s1q, s1c, s1d, s1f, s1t, s1x, nf, na, np, eub);

			// Initialize the inverse Hessian matrix to identity matrix.
			bfh[o0 = gid] = 1.0f;
			for (j = 1; j < nv; ++j)
			{
				for (i = 0; i < j; ++i)
				{
					bfh[o0 += gds] = 0.0f;
				}
				bfh[o0 += gds] = 1.0f;
			}
		}
		else
		{
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
	}

	// Accept x1 according to Metropolis criteria.
	if (s1e[gid] < s0e[gid])
	{
		o0 = gid;
		s0e[o0] = s1e[o0];
//		for (i = 1; i < nv + 2; ++i)
		for (i = -1 - nv; i < 0; ++i)
		{
			o0 += gds;
			s0e[o0] = s1e[o0];
		}
	}
}

int cu_mc_kernel::initialize(const int tid, const vector<float>& h_sf_e, const vector<float>& h_sf_d, const size_t h_sf_ns, const float* h_corner0, const float* h_corner1, const int* h_num_probes, const float h_granularity_inverse, const int num_mc_tasks, const int h_ng, const unsigned long h_seed)
{
	cudaSetDevice(device_id);
	this->num_mc_tasks = num_mc_tasks;

	// Initialize scoring function.
	const size_t sf_bytes = sizeof(float) * h_sf_e.size();
	checkCudaErrors(cudaMalloc(&d_sf_e, sf_bytes));
	checkCudaErrors(cudaMalloc(&d_sf_d, sf_bytes));
	checkCudaErrors(cudaMemcpy(d_sf_e, &h_sf_e.front(), sf_bytes, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_sf_d, &h_sf_d.front(), sf_bytes, cudaMemcpyHostToDevice));
	assert(sizeof(c_sf_e)  == sizeof(d_sf_e));
	assert(sizeof(c_sf_d)  == sizeof(d_sf_d));
//	assert(sizeof(c_sf_ns) == sizeof(h_sf_ns));
	checkCudaErrors(cudaMemcpyToSymbol(c_sf_e,  &d_sf_e,  sizeof(c_sf_e )));
	checkCudaErrors(cudaMemcpyToSymbol(c_sf_d,  &d_sf_d,  sizeof(c_sf_d )));
	checkCudaErrors(cudaMemcpyToSymbol(c_sf_ns, &h_sf_ns, sizeof(c_sf_ns)));

	// Initialize receptor.
	assert(sizeof(c_corner0) == sizeof(float) * 3);
	assert(sizeof(c_corner1) == sizeof(float) * 3);
	assert(sizeof(c_num_probes) == sizeof(int) * 3);
	assert(sizeof(c_granularity_inverse) == sizeof(h_granularity_inverse));
	assert(sizeof(c_ng) == sizeof(h_ng));
	checkCudaErrors(cudaMemcpyToSymbol(c_corner0, h_corner0, sizeof(c_corner0)));
	checkCudaErrors(cudaMemcpyToSymbol(c_corner1, h_corner1, sizeof(c_corner1)));
	checkCudaErrors(cudaMemcpyToSymbol(c_num_probes, h_num_probes, sizeof(c_num_probes)));
	checkCudaErrors(cudaMemcpyToSymbol(c_granularity_inverse, &h_granularity_inverse, sizeof(c_granularity_inverse)));
	checkCudaErrors(cudaMemcpyToSymbol(c_ng, &h_ng, sizeof(c_ng)));
	assert(sizeof(d_maps) == sizeof(float*) * sf_n);
	memset(d_maps, 0, sizeof(d_maps));

	// Initialize seed.
	assert(sizeof(c_seed) == sizeof(h_seed));
	checkCudaErrors(cudaMemcpyToSymbol(c_seed, &h_seed, sizeof(c_seed)));
	return 0;
}

int cu_mc_kernel::update(const int tid, const vector<vector<float> > h_maps, const vector<size_t>& xs)
{
	cudaSetDevice(device_id);
	const size_t map_bytes = sizeof(float) * h_maps[xs.front()].size();
	for (int i = 0; i < xs.size(); ++i)
	{
		const size_t t = xs[i];
		float* d_m;
		checkCudaErrors(cudaMalloc(&d_m, map_bytes));
		checkCudaErrors(cudaMemcpy(d_m, &h_maps[t].front(), map_bytes, cudaMemcpyHostToDevice));
		d_maps[t] = d_m;
	}
	assert(sizeof(c_maps) == sizeof(d_maps));
	checkCudaErrors(cudaMemcpyToSymbol(c_maps, d_maps, sizeof(c_maps)));
	return 0;
}

void cu_mc_kernel::launch(vector<float>& h_ex, const vector<int>& h_lig, const int nv, const int nf, const int na, const int np)
{
	cudaSetDevice(device_id);
	// Copy ligand content from host memory to device memory.
	const size_t lig_bytes = sizeof(int) * h_lig.size();
	int* d_lig;
	checkCudaErrors(cudaMalloc(&d_lig, lig_bytes));
	checkCudaErrors(cudaMemcpy(d_lig, &h_lig.front(), lig_bytes, cudaMemcpyHostToDevice));

	// Allocate device memory for variables. 3 * (nt + 1) is sufficient for t because the torques of inactive frames are always zero.
	const size_t var_bytes = sizeof(float) * ((1 + nv + 1 + nv + 3 * nf + 4 * nf + 3 * na + 3 * na + 3 * nf + 3 * nf) * 3 + (nv * (nv + 1) >> 1) + nv * 3) * num_mc_tasks;
	float* d_s0;
	checkCudaErrors(cudaMalloc(&d_s0, var_bytes));
	checkCudaErrors(cudaMemset(d_s0, 0, var_bytes));

	// Invoke CUDA kernel.
	mc<<<(num_mc_tasks - 1) / 32 + 1, 32, lig_bytes>>>(d_s0, d_lig, nv, nf, na, np);
	checkCudaErrors(cudaGetLastError());
	checkCudaErrors(cudaDeviceSynchronize());

	// Copy e and x from device memory to host memory.
	const size_t ex_size = (1 + nv + 1) * num_mc_tasks;
	const size_t ex_bytes = sizeof(float) * ex_size;
	h_ex.resize(ex_size);
	checkCudaErrors(cudaMemcpy(&h_ex.front(), d_s0, ex_bytes, cudaMemcpyDeviceToHost));

	// Free device memory.
	checkCudaErrors(cudaFree(d_s0));
	checkCudaErrors(cudaFree(d_lig));
}

cu_mc_kernel::~cu_mc_kernel()
{
	cudaSetDevice(device_id);
	for (size_t t = 0; t < sf_n; ++t)
	{
		float* const d_m = d_maps[t];
		if (d_m) checkCudaErrors(cudaFree(d_m));
	}
	checkCudaErrors(cudaFree(d_sf_d));
	checkCudaErrors(cudaFree(d_sf_e));
	checkCudaErrors(cudaDeviceReset());
}
