#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include "cu_helper.h"
#include "io_service_pool.hpp"
#include "receptor.hpp"
#include "ligand.hpp"
#include "utility.hpp"
#include "random_forest.hpp"
#include "log.hpp"
using namespace std;
using namespace boost::filesystem;

class safe_function
{
public:
	void operator()(function<void(void)>&& f)
	{
		lock_guard<mutex> guard(m);
		f();
	}
private:
	mutex m;
};

template <typename T>
class safe_counter
{
public:
	void init(const T z)
	{
		n = z;
		i = 0;
	}
	void increment()
	{
		lock_guard<mutex> guard(m);
		if (++i == n) cv.notify_one();
	}
	void wait()
	{
		unique_lock<mutex> lock(m);
		if (i < n) cv.wait(lock);
	}
private:
	mutex m;
	condition_variable cv;
	T n;
	T i;
};

template <typename T>
class safe_vector : public vector<T>
{
public:
	using vector<T>::vector;
	void safe_push_back(const T x)
	{
		lock_guard<mutex> guard(m);
		this->push_back(x);
		cv.notify_one();
	}
	T safe_pop_back()
	{
		unique_lock<mutex> lock(m);
		if (this->empty()) cv.wait(lock);
		const T x = this->back();
		this->pop_back();
		return x;
	}
private:
	mutex m;
	condition_variable cv;
};

int main(int argc, char* argv[])
{
	path receptor_path, input_folder_path, output_folder_path, log_path;
	array<float, 3> center, size;
	size_t seed, num_threads, num_trees, num_mc_tasks, num_bfgs_iterations, max_conformations;
	float granularity;

	// Parse program options in a try/catch block.
	try
	{
		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const path default_log_path = "log.csv";
		const size_t default_seed = chrono::system_clock::now().time_since_epoch().count();
		const size_t default_num_threads = thread::hardware_concurrency();
		const size_t default_num_trees = 128;
		const size_t default_num_mc_tasks = 256;
		const size_t default_num_bfgs_iterations = 300;
		const size_t default_max_conformations = 9;
		const  float default_granularity = 0.15625f;

		// Set up options description.
		using namespace boost::program_options;
		options_description input_options("input (required)");
		input_options.add_options()
			("receptor", value<path>(&receptor_path)->required(), "receptor in PDBQT format")
			("input_folder", value<path>(&input_folder_path)->required(), "folder of ligands in PDBQT format")
			("center_x", value<float>(&center[0])->required(), "x coordinate of the search space center")
			("center_y", value<float>(&center[1])->required(), "y coordinate of the search space center")
			("center_z", value<float>(&center[2])->required(), "z coordinate of the search space center")
			("size_x", value<float>(&size[0])->required(), "size in the x dimension in Angstrom")
			("size_y", value<float>(&size[1])->required(), "size in the y dimension in Angstrom")
			("size_z", value<float>(&size[2])->required(), "size in the z dimension in Angstrom")
			;
		options_description output_options("output (optional)");
		output_options.add_options()
			("output_folder", value<path>(&output_folder_path)->default_value(default_output_folder_path), "folder of output models in PDBQT format")
			("log", value<path>(&log_path)->default_value(default_log_path), "log file")
			;
		options_description miscellaneous_options("options (optional)");
		miscellaneous_options.add_options()
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "number of worker threads to use")
			("trees", value<size_t>(&num_trees)->default_value(default_num_trees), "number of trees in random forest")
			("tasks", value<size_t>(&num_mc_tasks)->default_value(default_num_mc_tasks), "number of Monte Carlo tasks for global search")
			("generations", value<size_t>(&num_bfgs_iterations)->default_value(default_num_bfgs_iterations), "number of generations in BFGS")
			("max_conformations", value<size_t>(&max_conformations)->default_value(default_max_conformations), "number of binding conformations to write")
			("granularity", value<float>(&granularity)->default_value(default_granularity), "density of probe atoms of grid maps")
			("help", "help information")
			("version", "version information")
			("config", value<path>(), "options can be loaded from a configuration file")
			;
		options_description all_options;
		all_options.add(input_options).add(output_options).add(miscellaneous_options);

		// Parse command line arguments.
		variables_map vm;
		store(parse_command_line(argc, argv, all_options), vm);

		// If no command line argument is supplied or help is requested, print the usage and exit.
		if (argc == 1 || vm.count("help"))
		{
			cout << all_options;
			return 0;
		}

		// If version is requested, print the version and exit.
		if (vm.count("version"))
		{
			cout << "3.0.0" << endl;
			return 0;
		}

		// If a configuration file is presented, parse it.
		if (vm.count("config"))
		{
			boost::filesystem::ifstream config_file(vm["config"].as<path>());
			store(parse_config_file(config_file, all_options), vm);
		}

		// Notify the user of parsing errors, if any.
		vm.notify();

		// Validate receptor.
		if (!is_regular_file(receptor_path))
		{
			cerr << "Receptor " << receptor_path << " does not exist or is not a regular file" << endl;
			return 1;
		}

		// Validate input_folder.
		if (!is_directory(input_folder_path))
		{
			cerr << "Input folder " << input_folder_path << " does not exist or is not a directory" << endl;
			return 1;
		}

		// Validate output_folder.
		if (exists(output_folder_path))
		{
			if (!is_directory(output_folder_path))
			{
				cerr << "Output folder " << output_folder_path << " is not a directory" << endl;
				return 1;
			}
		}
		else
		{
			if (!create_directories(output_folder_path))
			{
				cerr << "Failed to create output folder " << output_folder_path << endl;
				return 1;
			}
		}
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}

	cout << "Creating an io service pool of " << num_threads << " worker threads for host" << endl;
	io_service_pool ioh(num_threads);
	safe_counter<size_t> cnt;

	cout << "Precalculating a scoring function of " << scoring_function::n << " atom types in parallel" << endl;
	scoring_function sf;
	cnt.init(sf.n * (sf.n + 1) >> 1);
	for (size_t t2 = 0; t2 < sf.n; ++t2)
	for (size_t t1 = 0; t1 <=  t2; ++t1)
	{
		ioh.post([&,t1,t2]()
		{
			sf.precalculate(t1, t2);
			cnt.increment();
		});
	}
	cnt.wait();

	cout << "Parsing receptor " << receptor_path << endl;
	receptor rec(receptor_path, center, size, granularity);

	cout << "Detecting CUDA devices and just-in-time compiling modules" << endl;
	checkCudaErrors(cuInit(0));
	int num_devices;
	checkCudaErrors(cuDeviceGetCount(&num_devices));
	if (!num_devices)
	{
		cerr << "No CUDA devices detected" << endl;
		return 2;
	}
	cout << "DV               Name  CC SM GMEM(MB) SMEM(KB) CMEM(KB) MAPHOST ECC TIMEOUT MODE" << endl;
	vector<int> can_map_host_memory(num_devices);
	vector<CUcontext> contexts(num_devices);
	vector<CUfunction> functions(num_devices);
	vector<CUdeviceptr> mps(num_devices);
	for (int dev = 0; dev < num_devices; ++dev)
	{
		// Get a device handle from an ordinal.
		CUdevice device;
		checkCudaErrors(cuDeviceGet(&device, dev));

		// Get and print device attributes.
		char name[256];
		size_t totalGlobalMem;
		int major;
		int minor;
		int multiProcessorCount;
		int sharedMemPerBlock;
		int totalConstMem;
		int canMapHostMemory;
		int ECCEnabled;
		int kernelExecTimeoutEnabled;
		int computeMode;
		checkCudaErrors(cuDeviceGetName(name, sizeof(name), device));
		checkCudaErrors(cuDeviceTotalMem(&totalGlobalMem, device));
		checkCudaErrors(cuDeviceGetAttribute(&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, device));
		checkCudaErrors(cuDeviceGetAttribute(&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, device));
		checkCudaErrors(cuDeviceGetAttribute(&multiProcessorCount, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, device));
		checkCudaErrors(cuDeviceGetAttribute(&sharedMemPerBlock, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, device));
		checkCudaErrors(cuDeviceGetAttribute(&totalConstMem, CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY, device));
		checkCudaErrors(cuDeviceGetAttribute(&canMapHostMemory, CU_DEVICE_ATTRIBUTE_CAN_MAP_HOST_MEMORY, device));
		checkCudaErrors(cuDeviceGetAttribute(&ECCEnabled, CU_DEVICE_ATTRIBUTE_ECC_ENABLED, device));
		checkCudaErrors(cuDeviceGetAttribute(&kernelExecTimeoutEnabled, CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT, device));
		checkCudaErrors(cuDeviceGetAttribute(&computeMode, CU_DEVICE_ATTRIBUTE_COMPUTE_MODE, device));
		cout << setw(2) << dev << setw(19) << name << setw(2) << major << '.' << minor << setw(3) << multiProcessorCount << setw(9) << totalGlobalMem / 1048576 << setw(9) << sharedMemPerBlock / 1024 << setw(9) << totalConstMem / 1024 << setw(8) << canMapHostMemory << setw(4) << ECCEnabled << setw(8) << kernelExecTimeoutEnabled << setw(5) << computeMode << endl;

		// Save the device attribute of host memory mapping capability.
		can_map_host_memory[dev] = canMapHostMemory;

		// Create a context for the current device.
		checkCudaErrors(cuCtxCreate(&contexts[dev], CU_CTX_SCHED_AUTO/*CU_CTX_SCHED_YIELD*/ | (canMapHostMemory ? CU_CTX_MAP_HOST : 0), device));
//		checkCudaErrors(cuCtxSetCacheConfig(CU_FUNC_CACHE_PREFER_L1));
//		checkCudaErrors(cuCtxSetSharedMemConfig(CU_SHARED_MEM_CONFIG_EIGHT_BYTE_BANK_SIZE));

		// Initialize just-in-time compilation options.
//		const unsigned int numOptions = 6;
//		CUjit_option* options = (CUjit_option*)malloc(sizeof(CUjit_option) * numOptions);
//		void** optionValues = (void**)malloc(sizeof(void*) * numOptions);
//		options[0] = CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES;
//		size_t/*unsigned int*/ info_log_buffer_size_bytes = 1024;
//		optionValues[0] = (void*)info_log_buffer_size_bytes;
//		options[1] = CU_JIT_INFO_LOG_BUFFER;
//		char* info_log_buffer = (char*)malloc(sizeof(char) * info_log_buffer_size_bytes);
//		optionValues[1] = info_log_buffer;
//		options[2] = CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES;
//		size_t/*unsigned int*/ error_log_buffer_size_bytes = 1024;
//		optionValues[2] = (void*)error_log_buffer_size_bytes;
//		options[3] = CU_JIT_ERROR_LOG_BUFFER;
//		char* error_log_buffer = (char*)malloc(sizeof(char) * error_log_buffer_size_bytes);
//		optionValues[3] = error_log_buffer;
//		options[4] = CU_JIT_LOG_VERBOSE;
//		size_t/*int*/ log_verbose = 1;
//		optionValues[4] = (void*)log_verbose;
//		options[5] = CU_JIT_MAX_REGISTERS;
//		size_t/*unsigned int*/ max_registers = 32;
//		optionValues[5] = (void*)max_registers;
//		options[6] = CU_JIT_CACHE_MODE;
//		CUjit_cacheMode cache_mode = CU_JIT_CACHE_OPTION_NONE; // CU_JIT_CACHE_OPTION_CG, CU_JIT_CACHE_OPTION_CA
//		optionValues[6] = (void*)cache_mode;
//		options[7] = CU_JIT_GENERATE_DEBUG_INFO;
//		int generated_debug_info = 1;
//		optionValues[7] = (void*)generated_debug_info;
//		options[8] = CU_JIT_GENERATE_LINE_INFO;
//		int generated_line_info = 1;
//		optionValues[8] = (void*)generated_line_info;

		// Load the Monte Carlo module into the current context.
		CUmodule module;
		checkCudaErrors(cuModuleLoad(&module, "/home/hjli/idock/bin/idock.fatbin"));
//		checkCudaErrors(cuModuleLoadDataEx(&module, source, numOptions, options, optionValues));
//		printf("%s\n", info_log_buffer);
//		printf("%s\n", error_log_buffer);
//		free(error_log_buffer);
//		free(info_log_buffer);
//		free(optionValues);
//		free(options);

		// Get functions from module.
		checkCudaErrors(cuModuleGetFunction(&functions[dev], module, "monte_carlo"));

		// Get symbols from module.
		CUdeviceptr sfec;
		CUdeviceptr sfdc;
		CUdeviceptr sfsc;
		CUdeviceptr cr0c;
		CUdeviceptr cr1c;
		CUdeviceptr nprc;
		CUdeviceptr gric;
		CUdeviceptr mpsc;
		CUdeviceptr nbic;
		CUdeviceptr sedc;
		size_t sfes;
		size_t sfds;
		size_t sfss;
		size_t cr0s;
		size_t cr1s;
		size_t nprs;
		size_t gris;
		size_t mpss;
		size_t nbis;
		size_t seds;
		checkCudaErrors(cuModuleGetGlobal(&sfec, &sfes, module, "sfe")); //   8 const float*
		checkCudaErrors(cuModuleGetGlobal(&sfdc, &sfds, module, "sfd")); //   8 const float*
		checkCudaErrors(cuModuleGetGlobal(&sfsc, &sfss, module, "sfs")); //   4 int
		checkCudaErrors(cuModuleGetGlobal(&cr0c, &cr0s, module, "cr0")); //  12 float3
		checkCudaErrors(cuModuleGetGlobal(&cr1c, &cr1s, module, "cr1")); //  12 float3
		checkCudaErrors(cuModuleGetGlobal(&nprc, &nprs, module, "npr")); //  12 int3
		checkCudaErrors(cuModuleGetGlobal(&gric, &gris, module, "gri")); //   4 float
		checkCudaErrors(cuModuleGetGlobal(&mpsc, &mpss, module, "mps")); // 120 conat float* [15]
		checkCudaErrors(cuModuleGetGlobal(&nbic, &nbis, module, "nbi")); //   4 int
		checkCudaErrors(cuModuleGetGlobal(&sedc, &seds, module, "sed")); //   8 unsigned long

		// Initialize symbols for scoring function.
		CUdeviceptr sfed;
		CUdeviceptr sfdd;
		const int sfsh = sf.ns;
		assert(sfes == sizeof(sfed));
		assert(sfds == sizeof(sfdd));
		assert(sfss == sizeof(sfsh));
		const size_t sfe_bytes = sizeof(float) * sf.e.size();
		const size_t sfd_bytes = sizeof(float) * sf.d.size();
		checkCudaErrors(cuMemAlloc(&sfed, sfe_bytes));
		checkCudaErrors(cuMemAlloc(&sfdd, sfd_bytes));
		checkCudaErrors(cuMemcpyHtoD(sfed, sf.e.data(), sfe_bytes));
		checkCudaErrors(cuMemcpyHtoD(sfdd, sf.d.data(), sfd_bytes));
		checkCudaErrors(cuMemcpyHtoD(sfec, &sfed, sfes));
		checkCudaErrors(cuMemcpyHtoD(sfdc, &sfdd, sfds));
		checkCudaErrors(cuMemcpyHtoD(sfsc, &sfsh, sfss));

		// Initialize symbols for receptor.
		assert(cr0s == sizeof(rec.corner0));
		assert(cr1s == sizeof(rec.corner1));
		assert(nprs == sizeof(rec.num_probes));
		assert(gris == sizeof(rec.granularity_inverse));
		assert(mpss == sizeof(float*) * sf.n);
		checkCudaErrors(cuMemcpyHtoD(cr0c, rec.corner0.data(), cr0s));
		checkCudaErrors(cuMemcpyHtoD(cr1c, rec.corner1.data(), cr1s));
		checkCudaErrors(cuMemcpyHtoD(nprc, rec.num_probes.data(), nprs));
		checkCudaErrors(cuMemcpyHtoD(gric, &rec.granularity_inverse, gris));
		mps[dev] = mpsc;

		// Initialize symbols for program control.
		const int nbih = num_bfgs_iterations;
		assert(nbis == sizeof(nbih));
		assert(seds == sizeof(seed));
		checkCudaErrors(cuMemcpyHtoD(nbic, &nbih, nbis));
		checkCudaErrors(cuMemcpyHtoD(sedc, &seed, seds));

		// Pop the current context.
		checkCudaErrors(cuCtxPopCurrent(NULL));
	}
	sf.clear();

	// Initialize a vector of idle devices.
	safe_vector<int> idle(num_devices);
	iota(idle.begin(), idle.end(), 0);

	cout << "Training a random forest of " << num_trees << " trees with seed " << seed << " in parallel" << endl;
	forest f(num_trees, seed);
	cnt.init(num_trees);
	for (size_t i = 0; i < num_trees; ++i)
	{
		ioh.post([&,i]()
		{
			f[i].train(5, f.u01_s);
			cnt.increment();
		});
	}
	cnt.wait();
	f.clear();

	cout << "Creating an io service pool of " << num_devices << " worker threads for device" << endl;
	io_service_pool iod(num_devices);
	safe_function safe_print;

	// Perform docking for each ligand in the input folder.
	log_engine log;
	cout.setf(ios::fixed, ios::floatfield);
	cout << "Executing " << num_mc_tasks << " optimization runs of " << num_bfgs_iterations << " BFGS iterations in parallel" << endl;
	cout << "   Index       Ligand Dv   pKd 1     2     3     4     5     6     7     8     9" << endl << setprecision(2);
	for (directory_iterator dir_iter(input_folder_path), const_dir_iter; dir_iter != const_dir_iter; ++dir_iter)
	{
		// Parse the ligand. Don't declare it const as it will be moved to the io service pool for device.
		ligand lig(dir_iter->path());

		// Find atom types that are presented in the current ligand but not presented in the grid maps.
		vector<size_t> xs;
		for (const atom& a : lig.atoms)
		{
			const size_t t = a.xs;
			if (rec.maps[t].empty() && find(xs.cbegin(), xs.cend(), t) == xs.cend())
			{
				xs.push_back(t);
				rec.maps[t].resize(rec.num_probes_product);
			}
		}

		// Create grid maps on the fly if necessary.
		if (xs.size())
		{
			// Precalculate p_offset.
			for (size_t t1 = 0; t1 < sf.n; ++t1)
			{
				vector<size_t>& p = rec.p_offset[t1];
				p.resize(xs.size());
				for (size_t i = 0; i < xs.size(); ++i)
				{
					const size_t t2 = xs[i];
					p[i] = sf.nr * mp(t1, t2);
				}
			}

			// Create grid maps in parallel.
			cnt.init(rec.num_probes[2]);
			for (size_t z = 0; z < rec.num_probes[2]; ++z)
			{
				ioh.post([&,z]()
				{
					rec.populate(sf, xs, z);
					cnt.increment();
				});
			}
			cnt.wait();

			// Copy grid maps from host memory to device memory.
			for (int dev = 0; dev < num_devices; ++dev)
			{
				checkCudaErrors(cuCtxPushCurrent(contexts[dev]));
				const size_t map_bytes = sizeof(float) * rec.maps[xs.front()].size();
				for (const auto t : xs)
				{
					CUdeviceptr mapd;
					checkCudaErrors(cuMemAlloc(&mapd, map_bytes));
					checkCudaErrors(cuMemcpyHtoD(mapd, rec.maps[t].data(), map_bytes));
					checkCudaErrors(cuMemcpyHtoD(mps[dev] + sizeof(CUdeviceptr) * t, &mapd, sizeof(mapd)));
				}
				checkCudaErrors(cuCtxPopCurrent(NULL));
			}
			cnt.wait();
		}

		// Wait until a device is ready for execution.
		const int dev = idle.safe_pop_back();

		// Move the ligand from main thread to the io service pool for device.
		iod.post(bind<void>([&,dev](ligand& lig)
		{
			// Push the context of the chosen device.
			checkCudaErrors(cuCtxPushCurrent(contexts[dev]));

			// Map or copy ligand from host memory to device memory.
			const size_t lig_bytes = sizeof(int) * (11 * lig.nf + lig.nf - 1 + 4 * lig.na + 3 * lig.np);
			int* ligh;
			checkCudaErrors(cuMemHostAlloc((void**)&ligh, lig_bytes, can_map_host_memory[dev] ? CU_MEMHOSTALLOC_DEVICEMAP : 0));
			lig.encode(ligh);
			CUdeviceptr ligd;
			if (can_map_host_memory[dev])
			{
				checkCudaErrors(cuMemHostGetDevicePointer(&ligd, ligh, 0));
			}
			else
			{
				checkCudaErrors(cuMemAlloc(&ligd, lig_bytes));
				checkCudaErrors(cuMemcpyHtoDAsync(ligd, ligh, lig_bytes, NULL));
			}

			// Allocate device memory for solutions.
			// 3 * (nt + 1) is sufficient for t because the torques of inactive frames are always zero.
			const size_t sln_elems = ((1 + lig.nv + 1 + lig.nv + 3 * lig.nf + 4 * lig.nf + 3 * lig.na + 3 * lig.na + 3 * lig.nf + 3 * lig.nf) * 3 + (lig.nv * (lig.nv + 1) >> 1) + lig.nv * 3) * num_mc_tasks;
			const size_t sln_bytes = sizeof(float) * sln_elems;
			CUdeviceptr slnd;
			checkCudaErrors(cuMemAlloc(&slnd, sln_bytes));
			checkCudaErrors(cuMemsetD32Async(slnd, 0, sln_elems, NULL));

			// Launch kernel.
			void* params[] = { &slnd, &ligd, &lig.nv, &lig.nf, &lig.na, &lig.np };
			checkCudaErrors(cuLaunchKernel(functions[dev], (num_mc_tasks - 1) / 32 + 1, 1, 1, 32, 1, 1, lig_bytes, NULL, params, NULL));

			// Copy conformations from device memory to host memory.
			float* cnfh;
			const size_t cnf_bytes = sizeof(float) * (1 + lig.nv + 1) * num_mc_tasks;
			checkCudaErrors(cuMemHostAlloc((void**)&cnfh, cnf_bytes, 0));
			checkCudaErrors(cuMemcpyDtoHAsync(cnfh, slnd, cnf_bytes, NULL));

			// Synchronize.
			checkCudaErrors(cuCtxSynchronize());

			// Free device memory.
			checkCudaErrors(cuMemFree(slnd));
			if (!can_map_host_memory[dev])
			{
				checkCudaErrors(cuMemFree(ligd));
			}

			// Write conformations.
			lig.write(cnfh, output_folder_path, max_conformations, num_mc_tasks, rec, f);

			// Output and save ligand stem and predicted affinities.
			safe_print([&]()
			{
				string stem = lig.p.stem().string();
				cout << setw(8) << log.size() + 1 << setw(13) << stem << setw(3) << dev << "  ";
				for_each(lig.affinities.cbegin(), lig.affinities.cbegin() + min<size_t>(lig.affinities.size(), 9), [](const float a)
				{
					cout << setw(6) << a;
				});
				cout << endl;
				log.push_back(new log_record(move(stem), move(lig.affinities)));
			});

			// Free host memory.
			checkCudaErrors(cuMemFreeHost(ligh));
			checkCudaErrors(cuMemFreeHost(cnfh));

			// Pop the context after use.
			checkCudaErrors(cuCtxPopCurrent(NULL));

			// Signal the main thread to post another task.
			idle.safe_push_back(dev);
		}, move(lig)));
	}

	// Wait until the io service pool for host has finished all its tasks.
	ioh.wait();

	// Wait until the io service pool for device has finished all its tasks.
	iod.wait();

	// Destroy contexts.
	for (auto& context : contexts)
	{
//		const float* const mpsc = mps[dev];
//		for (size_t t = 0; t < sf.n; ++t)
//		{
//			float* const mapd = mpsc[t];
//			if (mapd) checkCudaErrors(cuMemFree(mapd));
//		}
//		checkCudaErrors(cuMemFree(sfdd));
//		checkCudaErrors(cuMemFree(sfed));
		checkCudaErrors(cuCtxDestroy(context));
	}

	// Sort and write ligand log records to the log file.
	if (log.empty()) return 0;
	cout << "Writing log records of " << log.size() << " ligands to " << log_path << endl;
	log.sort();
	log.write(log_path);
}
