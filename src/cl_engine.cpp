#include <sstream>
#include <CL/cl.h>
#include "cl_helper.h"
#include "cl_engine.hpp"

cl_engine::cl_engine()
{
	cl_uint num_platforms;
	checkOclErrors(clGetPlatformIDs(0, NULL, &num_platforms));
	if (!num_platforms) return;
	vector<cl_platform_id> platforms(num_platforms);
	checkOclErrors(clGetPlatformIDs(num_platforms, platforms.data(), NULL));
	for (cl_uint k = 0; k < num_platforms; ++k)
	{
		const cl_platform_id platform = platforms[k];
		cl_uint num_devices;
		const cl_int error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
		if (error == CL_DEVICE_NOT_FOUND) continue; // Intel(R) OpenCL sets num_devices = 1 even if clGetDeviceIDs returns CL_DEVICE_NOT_FOUND.
		checkOclErrors(error);
		vector<cl_device_id> devices(num_devices);
		checkOclErrors(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, num_devices, devices.data(), NULL));
		for (cl_uint j = 0; j < num_devices; ++j)
		{
			const cl_device_id device = devices[j];
			char name[256];
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(name), name, NULL));
			char opencl_c_version[256];
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, sizeof(opencl_c_version), opencl_c_version, NULL));
			cl_uint max_compute_units;
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(max_compute_units), &max_compute_units, NULL));
			cl_ulong global_mem_size;
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(global_mem_size), &global_mem_size, NULL));
			cl_ulong local_mem_size;
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local_mem_size), &local_mem_size, NULL));
			cl_ulong max_constant_buffer_size;
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(max_constant_buffer_size), &max_constant_buffer_size, NULL));
			cl_bool error_correction_support;
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(error_correction_support), &error_correction_support, NULL));
			ostringstream oss;
			oss << name << ", " << opencl_c_version << ", " << max_compute_units << " CUs, " << global_mem_size / 1048576 << "MB GMEM, " << local_mem_size / 1024 << "KB LMEM, " << max_constant_buffer_size / 1024 << "KB CMEM, ECC " << (error_correction_support ? "ON" : "OFF");
			this->devices.push_back(oss.str());
		}
	}
}
