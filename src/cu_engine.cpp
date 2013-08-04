#include <sstream>
#include <cuda_runtime_api.h>
#include <helper_cuda.h>
#include "cu_engine.hpp"

cu_engine::cu_engine()
{
	int num_devices;
	checkCudaErrors(cudaGetDeviceCount(&num_devices));

	devices.reserve(num_devices);
	for (int i = 0; i < num_devices; ++i)
	{
        checkCudaErrors(cudaSetDevice(i));
        cudaDeviceProp deviceProp;
        checkCudaErrors(cudaGetDeviceProperties(&deviceProp, i));
		ostringstream oss;
		oss << deviceProp.name << ", SM " << deviceProp.major << '.' << deviceProp.minor << ", " << deviceProp.multiProcessorCount << " SMs, " << deviceProp.totalGlobalMem / 1048576 << "MB GMEM, " << deviceProp.sharedMemPerBlock / 1024 << "KB SMEM, " << deviceProp.totalConstMem / 1024 << "KB CMEM, ECC " << (deviceProp.ECCEnabled ? "ON" : "OFF") << ", TIMEOUT " << (deviceProp.kernelExecTimeoutEnabled ? "ON" : "OFF");
		devices.push_back(oss.str());
	}
}
