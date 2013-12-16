/*
The idock kernel for OpenCL uses the MWC64X random number generator.

MWC64X
http://cas.ee.ic.ac.uk/people/dt10/research/rngs-gpu-mwc64x.html
Copyright (c) 2011, David Thomas
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice,
	this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright
	notice, this list of conditions and the following disclaimer in the
	documentation and/or other materials provided with the distribution.
	* Neither the name of Imperial College London nor the names of its
	contributors may be used to endorse or promote products derived
	from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// Pre: a < M, b < M
// Post: r = (a + b) mod M
ulong AddMod64(ulong a, ulong b, ulong M)
{
	ulong v = a + b;
	if (v >= M || v < a) v -= M;
	return v;
}

// Pre: a < M, b < M
// Post: r = (a * b) mod M
ulong MulMod64(ulong a, ulong b, ulong M)
{
	ulong r = 0;
	while (a)
	{
		if (a & 1) r = AddMod64(r, b, M);
		b = AddMod64(b, b, M);
		a = a >> 1;
	}
	return r;
}

// Pre: a < M, e >= 0
// Post: r = (a ^ b) mod M
// This takes at most ~64^2 modular additions, so probably about 2^15 or so instructions on most architectures
ulong PowMod64(ulong a, ulong e, ulong M)
{
	ulong sqr = a, acc = 1;
	while (e)
	{
		if (e & 1) acc = MulMod64(acc, sqr, M);
		sqr = MulMod64(sqr, sqr, M);
		e = e >> 1;
	}
	return acc;
}

typedef struct { uint x; uint c; } mwc64x_state_t;

enum { A = 4294883355U };
enum { M = 18446383549859758079UL };
enum { B = 4077358422479273989UL };

void skip(mwc64x_state_t *s, ulong d)
{
	ulong m = PowMod64(A, d, M);
	ulong x = MulMod64(s->x * (ulong)A + s->c, m, M);
	s->x = x / A;
	s->c = x % A;
}

void seed(mwc64x_state_t *s, ulong baseOffset, ulong perStreamOffset)
{
	ulong d = baseOffset + get_global_id(0) * perStreamOffset;
	ulong m = PowMod64(A, d, M);
	ulong x = MulMod64(B, m, M);
	s->x = x / A;
	s->c = x % A;
}

uint next(mwc64x_state_t *s)
{
	uint X = s->x;
	uint C = s->c;
	uint r = X ^ C;
	uint Xn = A * X + C;
	uint carry = Xn < C;
	uint Cn = mad_hi(A, X, carry);
	s->x = Xn;
	s->c = Cn;
	return r;
}

// Avoid using Shared Local Memory on the Intel Xeon Phi coprocessor.
// Have at least 1000 WGs per NDRange to optimally utilize Phi.
// Use Array Notation with int32 Indices.

#define assert(arg)

inline
bool evaluate(__global float* e, __global float* g, __global float* a, __global float* q, __global float* c, __global float* d, __global float* f, __global float* t, __global const float* x, const int nf, const int na, const int np, const float eub, __local const int* shared, __global const float* sfe, __global const float* sfd, const int sfs, const float3 cr0, const float3 cr1, const int3 npr, const float gri, __global const float* const mps[15])
{
	const int gid = get_global_id(0);
	const int gds = get_global_size(0);
	const int gd3 = 3 * gds;
	const int gd4 = 4 * gds;

	__local const int* const act = shared;
	__local const int* const beg = &act[nf];
	__local const int* const end = &beg[nf];
	__local const int* const nbr = &end[nf];
	__local const int* const prn = &nbr[nf];
	__local const float* const yy0 = (__local const float*)&prn[nf];
	__local const float* const yy1 = &yy0[nf];
	__local const float* const yy2 = &yy1[nf];
	__local const float* const xy0 = &yy2[nf];
	__local const float* const xy1 = &xy0[nf];
	__local const float* const xy2 = &xy1[nf];
	__local const int* const brs = (__local const int*)&xy2[nf];
	__local const float* const co0 = (__local const float*)&brs[nf - 1];
	__local const float* const co1 = &co0[na];
	__local const float* const co2 = &co1[na];
	__local const int* const xst = (__local const int*)&co2[na];
	__local const int* const ip0 = &xst[na];
	__local const int* const ip1 = &ip0[np];
	__local const int* const ipp = &ip1[np];

	float y, y0, y1, y2, v0, v1, v2, c0, c1, c2, e000, e100, e010, e001, a0, a1, a2, ang, sng, r0, r1, r2, r3, vs, dr, f0, f1, f2, t0, t1, t2, d0, d1, d2;
	float q0, q1, q2, q3, q00, q01, q02, q03, q11, q12, q13, q22, q23, q33, m0, m1, m2, m3, m4, m5, m6, m7, m8;
	int i, j, k, b, w, i0, i1, i2, k0, k1, k2, z;
	__global const float* map;

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
			assert(fabs(q0*q0 + q1*q1 + q2*q2 + q3*q3 - 1.0f) < 2e-3f);
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
			if (i == beg[k])
			{
				c0 = y0;
				c1 = y1;
				c2 = y2;
			}
			else
			{
				// Calculate coordinate from transformation matrix and offset.
				v0 = co0[i];
				v1 = co1[i];
				v2 = co2[i];
				c0 = y0 + m0 * v0 + m1 * v1 + m2 * v2;
				c1 = y1 + m3 * v0 + m4 * v1 + m5 * v2;
				c2 = y2 + m6 * v0 + m7 * v1 + m8 * v2;

				// Store coordinate from registers into memory.
				c[i0] = c0;
				c[i1] = c1;
				c[i2] = c2;
			}

			// TODO: move conditional expression out to bypass short circuiting.
			// Penalize out-of-box case.
			if (c0 < cr0.x || cr1.x <= c0 || c1 < cr0.y || cr1.y <= c1 || c2 < cr0.z || cr1.z <= c2)
			{
				y += 10.0f;
				d[i0] = 0.0f;
				d[i1] = 0.0f;
				d[i2] = 0.0f;
				continue;
			}

			// Find the index of the current coordinate
			k0 = (int)((c0 - cr0.x) * gri);
			k1 = (int)((c1 - cr0.y) * gri);
			k2 = (int)((c2 - cr0.z) * gri);
			assert(k0 + 1 < npr.x);
			assert(k1 + 1 < npr.y);
			assert(k2 + 1 < npr.z);
			k0 = npr.x * (npr.y * k2 + k1) + k0;

			// Retrieve the grid map and lookup the value
			 map = mps[xst[i]];
			e000 = map[k0];
			e100 = map[k0 + 1];
			e010 = map[k0 + npr.x];
			e001 = map[k0 + npr.x * npr.y];
			y += e000;
			d[i0] = (e100 - e000) * gri;
			d[i1] = (e010 - e000) * gri;
			d[i2] = (e001 - e000) * gri;
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
			assert(fabs(a0*a0 + a1*a1 + a2*a2 - 1.0f) < 2e-3f);
			a[k0  = i * gd3 + gid] = a0;
			a[k0 += gds] = a1;
			a[k0 += gds] = a2;

			// Update q of BRANCH frame
			ang = x[w += gds] * 0.5f;
//			sng = sin(ang);
//			r0 = cos(ang);
			sng = sincos(ang, &r0);
			r1 = sng * a0;
			r2 = sng * a1;
			r3 = sng * a2;
			q00 = r0 * q0 - r1 * q1 - r2 * q2 - r3 * q3;
			q01 = r0 * q1 + r1 * q0 + r2 * q3 - r3 * q2;
			q02 = r0 * q2 - r1 * q3 + r2 * q0 + r3 * q1;
			q03 = r0 * q3 + r1 * q2 - r2 * q1 + r3 * q0;
			assert(fabs(q00*q00 + q01*q01 + q02*q02 + q03*q03 - 1.0f) < 2e-3f);
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
		if (vs < 64.0f)
		{
			j = ipp[i] + (int)(sfs * vs);
			y += sfe[j];
			dr = sfd[j];
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

__kernel //__attribute__((reqd_work_group_size(X, Y, Z))) // X <= 16 (i.e. half warp or quarter wavefront) informs the compiler to optimize out barrier. Compile-time work group size helps the compiler to optimize register allocation.
void monte_carlo(__global float* const restrict s0e, __global const int* const restrict lig, const int nv, const int nf, const int na, const int np, const int nbi, __local int* const shared, __global const float* const sfe, __global const float* const sfd, const int sfs, const float3 cr0, const float3 cr1, const int3 npr, const float gri, __global const float* const x00, __global const float* const x01, __global const float* const x02, __global const float* const x03, __global const float* const x04, __global const float* const x05, __global const float* const x06, __global const float* const x07, __global const float* const x08, __global const float* const x09, __global const float* const x10, __global const float* const x11, __global const float* const x12, __global const float* const x13, __global const float* const x14)
{
	const int gid = get_global_id(0);
	const int gds = get_global_size(0);
	const int nls = 5; // Number of line search trials for determining step size in BFGS
	const float eub = 40.0f * na; // A conformation will be droped if its free energy is not better than e_upper_bound.
	__global float* const s0x = &s0e[gds];
	__global float* const s0g = &s0x[(nv + 1) * gds];
	__global float* const s0a = &s0g[nv * gds];
	__global float* const s0q = &s0a[3 * nf * gds];
	__global float* const s0c = &s0q[4 * nf * gds];
	__global float* const s0d = &s0c[3 * na * gds];
	__global float* const s0f = &s0d[3 * na * gds];
	__global float* const s0t = &s0f[3 * nf * gds];
	__global float* const s1e = &s0t[3 * nf * gds];
	__global float* const s1x = &s1e[gds];
	__global float* const s1g = &s1x[(nv + 1) * gds];
	__global float* const s1a = &s1g[nv * gds];
	__global float* const s1q = &s1a[3 * nf * gds];
	__global float* const s1c = &s1q[4 * nf * gds];
	__global float* const s1d = &s1c[3 * na * gds];
	__global float* const s1f = &s1d[3 * na * gds];
	__global float* const s1t = &s1f[3 * nf * gds];
	__global float* const s2e = &s1t[3 * nf * gds];
	__global float* const s2x = &s2e[gds];
	__global float* const s2g = &s2x[(nv + 1) * gds];
	__global float* const s2a = &s2g[nv * gds];
	__global float* const s2q = &s2a[3 * nf * gds];
	__global float* const s2c = &s2q[4 * nf * gds];
	__global float* const s2d = &s2c[3 * na * gds];
	__global float* const s2f = &s2d[3 * na * gds];
	__global float* const s2t = &s2f[3 * nf * gds];
	__global float* const bfh = &s2t[3 * nf * gds];
	__global float* const bfp = &bfh[(nv*(nv+1)>>1) * gds];
	__global float* const bfy = &bfp[nv * gds];
	__global float* const bfm = &bfy[nv * gds];
	float rd0, rd1, rd2, rd3, rst;
	float sum, pg1, pga, pgc, alp, pg2, pr0, pr1, pr2, nrm, ang, sng, pq0, pq1, pq2, pq3, s1xq0, s1xq1, s1xq2, s1xq3, s2xq0, s2xq1, s2xq2, s2xq3, bpi;
	float yhy, yps, ryp, pco, bpj, bmj, ppj;
	int g, i, j, o0, o1, o2;
	mwc64x_state_t rng;
	__global const float* const mps[15] = { x00, x01, x02, x03, x04, x05, x06, x07, x08, x09, x10, x11, x12, x13, x14 };

#ifdef CL_LOCAL
	// Load ligand into local memory.
	g = 11 * nf + nf - 1 + 4 * na + 3 * np;
	o0 = get_local_id(0);
	for (i = 0, j = (g - 1) / get_local_size(0); i < j; ++i)
	{
		shared[o0] = lig[o0];
		o0 += get_local_size(0);
	}
	if (o0 < g)
	{
		shared[o0] = lig[o0];
	}
	barrier(CLK_LOCAL_MEM_FENCE);
#else
#endif

	// Randomize s0x.
	seed(&rng, 0, 9999);
	rd0 = next(&rng);
	s0x[o0  = gid] = rd0 * cr1.x + (1 - rd0) * cr0.x;
	rd0 = next(&rng);
	s0x[o0 += gds] = rd0 * cr1.y + (1 - rd0) * cr0.y;
	rd0 = next(&rng);
	s0x[o0 += gds] = rd0 * cr1.z + (1 - rd0) * cr0.z;
	rd0 = next(&rng);
	rd1 = next(&rng);
	rd2 = next(&rng);
	rd3 = next(&rng);
	rst = rsqrt(rd0*rd0 + rd1*rd1 + rd2*rd2 + rd3*rd3);
	s0x[o0 += gds] = rd0 * rst;
	s0x[o0 += gds] = rd1 * rst;
	s0x[o0 += gds] = rd2 * rst;
	s0x[o0 += gds] = rd3 * rst;
	for (i = 6; i < nv; ++i)
	{
		s0x[o0 += gds] = next(&rng);
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
	evaluate(s0e, s0g, s0a, s0q, s0c, s0d, s0f, s0t, s0x, nf, na, np, eub, shared, sfe, sfd, sfs, cr0, cr1, npr, gri, mps);

	// Repeat for a number of generations.
	for (g = 0; g < nbi; ++g)
	{
		// Mutate s0x into s1x
		o0  = gid;
		s1x[o0] = s0x[o0] + next(&rng);
		o0 += gds;
		s1x[o0] = s0x[o0] + next(&rng);
		o0 += gds;
		s1x[o0] = s0x[o0] + next(&rng);
//		for (i = 3; i < nv + 1; ++i)
		for (i = 2 - nv; i < 0; ++i)
		{
			o0 += gds;
			s1x[o0] = s0x[o0];
		}
		evaluate(s1e, s1g, s1a, s1q, s1c, s1d, s1f, s1t, s1x, nf, na, np, eub, shared, sfe, sfd, sfs, cr0, cr1, npr, gri, mps);

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
				assert(fabs(s1xq0*s1xq0 + s1xq1*s1xq1 + s1xq2*s1xq2 + s1xq3*s1xq3 - 1.0f) < 2e-3f);
				nrm = sqrt(pr0*pr0 + pr1*pr1 + pr2*pr2);
				ang = 0.5f * alp * nrm;
//				sng = sin(ang) / nrm;
//				pq0 = cos(ang);
				sng = sincos(ang, &pq0) / nrm;
				pq1 = sng * pr0;
				pq2 = sng * pr1;
				pq3 = sng * pr2;
				assert(fabs(pq0*pq0 + pq1*pq1 + pq2*pq2 + pq3*pq3 - 1.0f) < 2e-3f);
				s2xq0 = pq0 * s1xq0 - pq1 * s1xq1 - pq2 * s1xq2 - pq3 * s1xq3;
				s2xq1 = pq0 * s1xq1 + pq1 * s1xq0 + pq2 * s1xq3 - pq3 * s1xq2;
				s2xq2 = pq0 * s1xq2 - pq1 * s1xq3 + pq2 * s1xq0 + pq3 * s1xq1;
				s2xq3 = pq0 * s1xq3 + pq1 * s1xq2 - pq2 * s1xq1 + pq3 * s1xq0;
				assert(fabs(s2xq0*s2xq0 + s2xq1*s2xq1 + s2xq2*s2xq2 + s2xq3*s2xq3 - 1.0f) < 2e-3f);
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
				if (evaluate(s2e, s2g, s2a, s2q, s2c, s2d, s2f, s2t, s2x, nf, na, np, s1e[gid] + alp * pga, shared, sfe, sfd, sfs, cr0, cr1, npr, gri, mps))
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
