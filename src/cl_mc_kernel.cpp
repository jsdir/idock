#include <cassert>
#include <cstring>
#include "cl_helper.h"
#include "cl_mc_kernel.hpp"
using namespace std;

void cl_mc_kernel::initialize(const float* h_sf_e, const float* h_sf_d, const int h_sf_ns, const int h_sf_ne, const float* h_corner0, const float* h_corner1, const int* h_num_probes, const float h_granularity_inverse, const int num_mc_tasks, const int h_ng, const unsigned long h_seed)
{
	// Find an appropriate platform.
	cl_uint num_platforms;
	cl_platform_id* platforms;
	cl_platform_id platform;
	cl_device_id device;
	checkOclErrors(clGetPlatformIDs(0, NULL, &num_platforms));
//	if (num_platforms == 0)
//	{
//		cerr << "No OpenCL platform found" << endl;
//		return 1;
//	}
	platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id) * num_platforms);
	checkOclErrors(clGetPlatformIDs(num_platforms, platforms, NULL));
	for (cl_uint i = 0; i < num_platforms; ++i)
	{
		// 0 AMD Accelerated Parallel Processing
		// 1 NVIDIA CUDA
		// 2 Intel(R) OpenCL
		char buffer[1024];
		checkOclErrors(clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, sizeof(buffer), buffer, NULL));
		if (!strcmp(buffer, "NVIDIA CUDA"))
		{
			platform = platforms[i];
			break;
		}
	}
	free(platforms);

	// Find a GPU device.
	checkOclErrors(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL));
	cl_context_properties context_properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0};
	context = clCreateContext(context_properties, 1, &device, NULL, NULL, &err);
	checkOclErrors(err);
	command_queue = clCreateCommandQueue(context, device, 0, &err);
	checkOclErrors(err);

	// JIT OpenCL kernel source.
	FILE* pFileStream = fopen("cl_mc_kernel.cl", "rb");
	fseek(pFileStream, 0, SEEK_END);
	size_t szSourceLength = ftell(pFileStream);
	fseek(pFileStream, 0, SEEK_SET);
	char* cSourceString = (char *)malloc(szSourceLength);
	fread(cSourceString, szSourceLength, 1, pFileStream);
	fclose(pFileStream);
	program = clCreateProgramWithSource(context, 1, (const char **)&cSourceString, &szSourceLength, &err);
	checkOclErrors(err);
	free(cSourceString);
	checkOclErrors(clBuildProgram(program, 0, NULL, "-cl-fast-relaxed-math -cl-denorms-are-zero", NULL, NULL));
	size_t build_log_size;
	checkOclErrors(clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &build_log_size));
	char* build_log = (char *)malloc(build_log_size);
	checkOclErrors(clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, build_log_size, build_log, NULL));
//	cout << build << endl;
	kernel = clCreateKernel(program, "bfgs", &err);
	checkOclErrors(err);

	this->num_mc_tasks = num_mc_tasks;

	// Initialize scoring function.
	d_sf_e = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * h_sf_ne, NULL, &err);
	checkOclErrors(err);
	d_sf_d = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * h_sf_ne, NULL, &err);
	checkOclErrors(err);
	checkOclErrors(clEnqueueWriteBuffer(command_queue, d_sf_e, CL_FALSE, 0, sizeof(float) * h_sf_ne, h_sf_e, 0, NULL, NULL));
	checkOclErrors(clEnqueueWriteBuffer(command_queue, d_sf_d, CL_FALSE, 0, sizeof(float) * h_sf_ne, h_sf_d, 0, NULL, NULL));
//	assert(sizeof(c_sf_e)  == sizeof(d_sf_e));
//	assert(sizeof(c_sf_d)  == sizeof(d_sf_d));
//	assert(sizeof(c_sf_ns) == sizeof(h_sf_ns));
//	checkOclErrors(cudaMemcpyToSymbol(c_sf_e,  &d_sf_e,  sizeof(c_sf_e )));
//	checkOclErrors(cudaMemcpyToSymbol(c_sf_d,  &d_sf_d,  sizeof(c_sf_d )));
//	checkOclErrors(cudaMemcpyToSymbol(c_sf_ns, &h_sf_ns, sizeof(c_sf_ns)));
//	c_sf_e = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(c_sf_e ), &d_sf_e, &err);

	// Initialize receptor.
//	assert(sizeof(c_corner0) == sizeof(float) * 3);
//	assert(sizeof(c_corner1) == sizeof(float) * 3);
//	assert(sizeof(c_num_probes) == sizeof(int) * 3);
//	assert(sizeof(c_granularity_inverse) == sizeof(h_granularity_inverse));
//	assert(sizeof(c_ng) == sizeof(h_ng));
//	checkOclErrors(cudaMemcpyToSymbol(c_corner0, h_corner0, sizeof(c_corner0)));
//	checkOclErrors(cudaMemcpyToSymbol(c_corner1, h_corner1, sizeof(c_corner1)));
//	checkOclErrors(cudaMemcpyToSymbol(c_num_probes, h_num_probes, sizeof(c_num_probes)));
//	checkOclErrors(cudaMemcpyToSymbol(c_granularity_inverse, &h_granularity_inverse, sizeof(c_granularity_inverse)));
//	checkOclErrors(cudaMemcpyToSymbol(c_ng, &h_ng, sizeof(c_ng)));
	assert(sizeof(d_maps) == sizeof(float*) * sf_n);
	memset(d_maps, 0, sizeof(d_maps));

	// Initialize seed.
//	assert(sizeof(c_seed) == sizeof(h_seed));
//	checkOclErrors(cudaMemcpyToSymbol(c_seed, &h_seed, sizeof(c_seed)));
}

void cl_mc_kernel::update(const vector<vector<float> > h_maps, const vector<size_t>& xs)
{
	const size_t map_bytes = sizeof(float) * h_maps[xs.front()].size();
	for (int i = 0; i < xs.size(); ++i)
	{
		const size_t t = xs[i];
		float* d_m;
//		checkOclErrors(cudaMalloc(&d_m, map_bytes));
//		checkOclErrors(cudaMemcpy(d_m, &h_maps[t].front(), map_bytes, cudaMemcpyHostToDevice));
//		d_maps[t] = d_m;
	}
//	assert(sizeof(c_maps) == sizeof(d_maps));
//	checkOclErrors(cudaMemcpyToSymbol(c_maps, d_maps, sizeof(c_maps)));
}

void cl_mc_kernel::launch(vector<float>& h_ex, const vector<int>& h_lig, const int nv, const int nf, const int na, const int np)
{
	// Copy ligand content from host memory to device memory.
	const size_t lig_bytes = sizeof(int) * h_lig.size();
	cl_mem d_lig = clCreateBuffer(context, CL_MEM_READ_ONLY/* | CL_MEM_COPY_HOST_PTR*/, lig_bytes, NULL/*&h_lig.front()*/, &err);
	checkOclErrors(err);
	checkOclErrors(clEnqueueWriteBuffer(command_queue, d_lig, CL_FALSE, 0, lig_bytes, &h_lig.front(), 0, NULL, NULL));

	// Allocate device memory for variables. 3 * (nt + 1) is sufficient for t because the torques of inactive frames are always zero.
	vector<float> h_s0(((1 + nv + 1 + nv + 3 * nf + 4 * nf + 3 * na + 3 * na + 3 * nf + 3 * nf) * 3 + (nv * (nv + 1) >> 1) + nv * 3) * num_mc_tasks, 0);
	const size_t s0_bytes = sizeof(float) * h_s0.size();
	cl_mem d_s0 = clCreateBuffer(context, CL_MEM_READ_WRITE, s0_bytes, NULL, &err);
	checkOclErrors(err);
	checkOclErrors(clEnqueueWriteBuffer(command_queue, d_s0, CL_FALSE, 0, s0_bytes, h_s0.data(), 0, NULL, NULL));

	// Invoke OpenCL kernel.
	checkOclErrors(clSetKernelArg(kernel, 0, sizeof(cl_mem), &d_s0));
	checkOclErrors(clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_lig));
	checkOclErrors(clSetKernelArg(kernel, 2, sizeof(int), &nv));
	checkOclErrors(clSetKernelArg(kernel, 3, sizeof(int), &nf));
	checkOclErrors(clSetKernelArg(kernel, 4, sizeof(int), &na));
	checkOclErrors(clSetKernelArg(kernel, 5, sizeof(int), &np));
	checkOclErrors(clSetKernelArg(kernel, 6, lig_bytes, NULL));
	const size_t lws = 32;
	const size_t gws = num_mc_tasks / lws;
//	checkOclErrors(clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &szMaxWorkgroupSize, NULL));
	checkOclErrors(clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &gws, &lws, 0, NULL, NULL));

	// Copy e and x from device memory to host memory.
	const size_t ex_size = (1 + nv + 1) * num_mc_tasks;
	const size_t ex_bytes = sizeof(float) * ex_size;
	h_ex.resize(ex_size);
	checkOclErrors(clEnqueueReadBuffer(command_queue, d_s0, CL_TRUE, 0, ex_bytes, &h_ex.front(), 0, NULL, NULL));

	// Free device memory.
	checkOclErrors(clReleaseMemObject(d_s0));
	checkOclErrors(clReleaseMemObject(d_lig));
}

cl_mc_kernel::~cl_mc_kernel()
{
	for (size_t t = 0; t < sf_n; ++t)
	{
		const cl_mem d_m = d_maps[t];
		if (d_m) checkOclErrors(clReleaseMemObject(d_m));
	}
	checkOclErrors(clReleaseMemObject(d_sf_d));
	checkOclErrors(clReleaseMemObject(d_sf_e));
	checkOclErrors(clReleaseKernel(kernel));
	checkOclErrors(clReleaseProgram(program));
	checkOclErrors(clReleaseCommandQueue(command_queue));
	checkOclErrors(clReleaseContext(context));
}
