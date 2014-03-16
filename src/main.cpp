#include <chrono>
#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include "io_service_pool.hpp"
#include "safe_counter.hpp"
#include "receptor.hpp"
#include "ligand.hpp"
#include "random_forest.hpp"
#include "log.hpp"

int main(int argc, char* argv[])
{
	path receptor_path, input_folder_path, output_folder_path, log_path;
	array<double, 3> center, size;
	size_t seed, num_threads, num_trees, num_mc_tasks, max_conformations;
	double grid_granularity;

	// Process program options.
	try
	{
		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const path default_log_path = "log.csv";
		const size_t default_seed = chrono::system_clock::now().time_since_epoch().count();
		const size_t default_num_threads = thread::hardware_concurrency();
		const size_t default_num_trees = 500;
		const size_t default_num_mc_tasks = 32;
		const size_t default_max_conformations = 9;
		const double default_grid_granularity = 0.15625;

		// Set up options description.
		using namespace boost::program_options;
		options_description input_options("input (required)");
		input_options.add_options()
			("receptor", value<path>(&receptor_path)->required(), "receptor in PDBQT format")
			("input_folder", value<path>(&input_folder_path)->required(), "folder of ligands in PDBQT format")
			("center_x", value<double>(&center[0])->required(), "x coordinate of the search space center")
			("center_y", value<double>(&center[1])->required(), "y coordinate of the search space center")
			("center_z", value<double>(&center[2])->required(), "z coordinate of the search space center")
			("size_x", value<double>(&size[0])->required(), "size in the x dimension in Angstrom")
			("size_y", value<double>(&size[1])->required(), "size in the y dimension in Angstrom")
			("size_z", value<double>(&size[2])->required(), "size in the z dimension in Angstrom")
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
			("max_conformations", value<size_t>(&max_conformations)->default_value(default_max_conformations), "maximum number of binding conformations to write")
			("granularity", value<double>(&grid_granularity)->default_value(default_grid_granularity), "density of probe atoms of grid maps")
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
			cout << "2.1" << endl;
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
		if (!exists(receptor_path))
		{
			cerr << "Receptor " << receptor_path << " does not exist" << endl;
			return 1;
		}
		if (!is_regular_file(receptor_path))
		{
			cerr << "Receptor " << receptor_path << " is not a regular file" << endl;
			return 1;
		}

		// Validate input_folder_path.
		if (!exists(input_folder_path))
		{
			cerr << "Input folder " << input_folder_path << " does not exist" << endl;
			return 1;
		}
		if (!is_directory(input_folder_path))
		{
			cerr << "Input folder " << input_folder_path << " is not a directory" << endl;
			return 1;
		}

		// Validate size.
		if (size[0] < receptor::Default_Partition_Granularity ||
			size[1] < receptor::Default_Partition_Granularity ||
			size[2] < receptor::Default_Partition_Granularity)
		{
			cerr << "Search space must be "
				 << receptor::Default_Partition_Granularity << "A x "
				 << receptor::Default_Partition_Granularity << "A x "
				 << receptor::Default_Partition_Granularity << "A or larger" << endl;
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

		// Validate log_path.
		if (is_directory(log_path))
		{
			cerr << "Option log " << log_path << " is a directory" << endl;
			return 1;
		}

		// Validate miscellaneous options.
		if (!num_threads)
		{
			cerr << "Option threads must be 1 or greater" << endl;
			return 1;
		}
		if (!num_mc_tasks)
		{
			cerr << "Option tasks must be 1 or greater" << endl;
			return 1;
		}
		if (!max_conformations)
		{
			cerr << "Option max_conformations must be 1 or greater" << endl;
			return 1;
		}
		if (grid_granularity <= 0)
		{
			cerr << "Option granularity must be positive" << endl;
			return 1;
		}
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}

	// Initialize a Mersenne Twister random number generator.
	cout << "Using random seed " << seed << endl;
	mt19937_64 rng(seed);

	// Initialize an io service pool and create worker threads for later use.
	cout << "Creating an io service pool of " << num_threads << " worker thread" << (num_threads == 1 ? "" : "s") << endl;
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;

	// Precalculate the scoring function in parallel.
	cout << "Precalculating scoring function in parallel" << endl;
	scoring_function sf;
	cnt.init(((sf.n + 1) * sf.n) >> 1);
	for (size_t t0 =  0; t0 < sf.n; ++t0)
	for (size_t t1 = t0; t1 < sf.n; ++t1)
	{
		io.post([&, t0, t1]()
		{
			sf.precalculate(t0, t1);
			cnt.increment();
		});
	}
	cnt.wait();
	sf.clear();

	// Parse the receptor.
	cout << "Parsing receptor " << receptor_path << endl;
	receptor rec(receptor_path, center, size, grid_granularity);
	const size_t num_gm_tasks = rec.num_probes[2];

	// Reserve storage for result containers. ptr_vector<T> is used for fast sorting.
	const size_t max_results = 20; // Maximum number of results obtained from a single Monte Carlo task.
	ptr_vector<ptr_vector<result>> result_containers;
	result_containers.resize(num_mc_tasks);
	for (size_t i = 0; i < num_mc_tasks; ++i)
	{
		result_containers[i].reserve(max_results);
	}
	ptr_vector<result> results;
	results.reserve(max_results * num_mc_tasks);

	cout << "Training a random forest of " << num_trees << " trees in parallel" << endl;
	forest f(num_trees, seed);
	cnt.init(num_trees);
	for (size_t i = 0; i < num_trees; ++i)
	{
		io.post([&, i]()
		{
			f[i].train(4, f.u01_s);
			cnt.increment();
		});
	}
	cnt.wait();
	f.clear();

	// Perform docking for each file in the ligand folder.
	cout << "Running " << num_mc_tasks << " Monte Carlo task" << (num_mc_tasks == 1 ? "" : "s") << " per ligand" << endl
	     << "   Index        Ligand Energy 1     2     3     4     5     6     7     8     9" << endl << setprecision(2);
	cout.setf(ios::fixed, ios::floatfield);
	log_engine log;
	for (directory_iterator dir_iter(input_folder_path), end_dir_iter; dir_iter != end_dir_iter; ++dir_iter)
	{
		// Filter files with .pdbqt extension name.
		const path input_ligand_path = dir_iter->path();
		if (input_ligand_path.extension() != ".pdbqt") continue;

		// Parse the ligand.
		ligand lig(input_ligand_path);

		// Find atom types that are presented in the current ligand but not presented in the grid maps.
		vector<size_t> xs;
		for (size_t t = 0; t < sf.n; ++t)
		{
			if (lig.xs[t] && rec.grid_maps[t].empty())
			{
				rec.grid_maps[t].resize(rec.num_probes_product);
				xs.push_back(t);
			}
		}

		// Create grid maps on the fly if necessary.
		if (xs.size())
		{
			// Populate the grid map task container.
			cnt.init(num_gm_tasks);
			for (size_t z = 0; z < num_gm_tasks; ++z)
			{
				io.post([&, z]()
				{
					rec.populate(xs, z, sf);
					cnt.increment();
				});
			}
			cnt.wait();
		}

		// Dump the ligand file stem.
		string stem = input_ligand_path.stem().string();
		cout << setw(8) << log.size() + 1 << setw(14) << stem << "   " << flush;

		// Populate the Monte Carlo task container.
		cnt.init(num_mc_tasks);
		for (size_t i = 0; i < num_mc_tasks; ++i)
		{
			assert(result_containers[i].empty());
			const size_t s = rng();
			io.post([&, i, s]()
			{
				lig.monte_carlo(result_containers[i], s, sf, rec);
				cnt.increment();
			});
		}
		cnt.wait();

		// Merge results from all the tasks into one single result container.
		assert(results.empty());
		const double required_square_error = 4.0 * lig.num_heavy_atoms; // Ligands with RMSD < 2.0 will be clustered into the same cluster.
		for (size_t i = 0; i < num_mc_tasks; ++i)
		{
			ptr_vector<result>& task_results = result_containers[i];
			const size_t num_task_results = task_results.size();
			for (size_t j = 0; j < num_task_results; ++j)
			{
				add_to_result_container(results, move(task_results[j]), required_square_error);
			}
			task_results.clear();
		}

		// If no conformation can be found, skip the current ligand and proceed with the next one.
		if (results.empty())
		{
			cout << endl;
			continue;
		}

		// Adjust free energy relative to the best conformation and flexibility.
		const size_t num_results = min<size_t>(results.size(), max_conformations);
		const result& best_result = results.front();
		const double best_result_intra_e = best_result.e - best_result.f;
		for (size_t i = 0; i < num_results; ++i)
		{
			results[i].e_nd = (results[i].e - best_result_intra_e) * lig.flexibility_penalty_factor;
		}

		// Write models to file.
		const path output_ligand_path = output_folder_path / input_ligand_path.filename();
		lig.write_models(output_ligand_path, results, num_results, rec, f, sf);

		// Display the free energies of the top 9 conformations.
		const size_t num_energies = min<size_t>(num_results, 9);
		for (size_t i = 0; i < num_energies; ++i)
		{
			cout << setw(6) << results[i].e_nd;
		}
		cout << endl;

		// Add a log record.
		vector<double> affinities(num_results);
		for (size_t i = 0; i < num_results; ++i)
		{
			affinities[i] = results[i].e_nd;
		}
		log.push_back(new log_record(move(stem), move(affinities)));

		// Clear the results of the current ligand.
		results.clear();
	}

	// Sort and write ligand log records to the log file.
	if (log.empty()) return 0;
	cout << "Writing log records of " << log.size() << " ligands to " << log_path << endl;
	log.sort();
	log.write(log_path);
}
