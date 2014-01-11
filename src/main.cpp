#include <chrono>
#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include "receptor.hpp"
#include "ligand.hpp"
#include "io_service_pool.hpp"
#include "grid_map_task.hpp"
#include "monte_carlo_task.hpp"
#include "utility.hpp"
//#include "random_forest.hpp"
#include "log.hpp"
using namespace std;

int main(int argc, char* argv[])
{
	path receptor_path, input_folder_path, output_folder_path, log_path;
	fl center_x, center_y, center_z, size_x, size_y, size_z;
	size_t num_threads, seed, num_mc_tasks, max_conformations;
	fl grid_granularity;

	// Process program options.
	try
	{
		using namespace boost::program_options;

		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const path default_log_path = "log.csv";
		const unsigned int concurrency = thread::hardware_concurrency();
		const size_t default_num_threads = concurrency;
		const size_t default_seed = chrono::system_clock::now().time_since_epoch().count();
		const size_t default_num_mc_tasks = 32;
		const size_t default_max_conformations = 9;
		const fl default_grid_granularity = 0.15625;

		options_description input_options("input (required)");
		input_options.add_options()
			("receptor", value<path>(&receptor_path)->required(), "receptor in PDBQT format")
			("input_folder", value<path>(&input_folder_path)->required(), "folder of ligands in PDBQT format")
			("center_x", value<fl>(&center_x)->required(), "x coordinate of the search space center")
			("center_y", value<fl>(&center_y)->required(), "y coordinate of the search space center")
			("center_z", value<fl>(&center_z)->required(), "z coordinate of the search space center")
			("size_x", value<fl>(&size_x)->required(), "size in the x dimension in Angstrom")
			("size_y", value<fl>(&size_y)->required(), "size in the y dimension in Angstrom")
			("size_z", value<fl>(&size_z)->required(), "size in the z dimension in Angstrom")
			;

		options_description output_options("output (optional)");
		output_options.add_options()
			("output_folder", value<path>(&output_folder_path)->default_value(default_output_folder_path), "folder of output models in PDBQT format")
			("log", value<path>(&log_path)->default_value(default_log_path), "log file")
			;

		options_description miscellaneous_options("options (optional)");
		miscellaneous_options.add_options()
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "number of worker threads to use")
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("tasks", value<size_t>(&num_mc_tasks)->default_value(default_num_mc_tasks), "number of Monte Carlo tasks for global search")
			("max_conformations", value<size_t>(&max_conformations)->default_value(default_max_conformations), "maximum number of binding conformations to write")
			("granularity", value<fl>(&grid_granularity)->default_value(default_grid_granularity), "density of probe atoms of grid maps")
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
			cerr << "Receptor " << receptor_path << " does not exist\n";
			return 1;
		}
		if (!is_regular_file(receptor_path))
		{
			cerr << "Receptor " << receptor_path << " is not a regular file\n";
			return 1;
		}

		// Validate input_folder_path.
		if (!exists(input_folder_path))
		{
			cerr << "Input folder " << input_folder_path << " does not exist\n";
			return 1;
		}
		if (!is_directory(input_folder_path))
		{
			cerr << "Input folder " << input_folder_path << " is not a directory\n";
			return 1;
		}

		// Validate size_x, size_y, size_z.
		if (size_x < box::Default_Partition_Granularity ||
		    size_y < box::Default_Partition_Granularity ||
		    size_z < box::Default_Partition_Granularity)
		{
			cerr << "Search space must be "
				 << box::Default_Partition_Granularity << "A x "
				 << box::Default_Partition_Granularity << "A x "
				 << box::Default_Partition_Granularity << "A or larger\n";
			return 1;
		}

		// Validate output_folder.
		if (exists(output_folder_path))
		{
			if (!is_directory(output_folder_path))
			{
				cerr << "Output folder " << output_folder_path << " is not a directory\n";
				return 1;
			}
		}
		else
		{
			if (!create_directories(output_folder_path))
			{
				cerr << "Failed to create output folder " << output_folder_path << '\n';
				return 1;
			}
		}

		// Validate log_path.
		if (is_directory(log_path))
		{
			cerr << "log path " << log_path << " is a directory\n";
			return 1;
		}

		// Validate miscellaneous options.
		if (!num_threads)
		{
			cerr << "Option threads must be 1 or greater\n";
			return 1;
		}
		if (!num_mc_tasks)
		{
			cerr << "Option tasks must be 1 or greater\n";
			return 1;
		}
		if (!max_conformations)
		{
			cerr << "Option max_conformations must be 1 or greater\n";
			return 1;
		}
		if (grid_granularity <= 0)
		{
			cerr << "Option granularity must be positive\n";
			return 1;
		}
	}
	catch (const exception& e)
	{
		cerr << e.what() << '\n';
		return 1;
	}

	// Initialize a Mersenne Twister random number generator.
	cout << "Using random seed " << seed << '\n';
	mt19937_64 eng(seed);

	// Initialize the search space of cuboid shape.
	const box b(vec3(center_x, center_y, center_z), vec3(size_x, size_y, size_z), grid_granularity);

	// Parse the receptor.
	cout << "Parsing receptor " << receptor_path << '\n';
	const receptor rec(receptor_path, b);

	// Reserve storage for task containers.
	const size_t num_gm_tasks = b.num_probes[0];

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

	// Initialize a vector of empty grid maps. Each grid map corresponds to an XScore atom type.
	vector<array3d<fl>> grid_maps(XS_TYPE_SIZE);
	vector<size_t> atom_types_to_populate;
	atom_types_to_populate.reserve(XS_TYPE_SIZE);

	// Initialize a thread pool and create worker threads for later use.
	cout << "Creating an io service pool of " << num_threads << " worker thread" << ((num_threads == 1) ? "" : "s") << '\n';
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;

	// Precalculate the scoring function in parallel.
	cout << "Precalculating scoring function in parallel ";
	scoring_function sf;
	{
		// Precalculate reciprocal square root values.
		vector<fl> rs(scoring_function::Num_Samples, 0);
		for (size_t i = 0; i < scoring_function::Num_Samples; ++i)
		{
			rs[i] = sqrt(i * scoring_function::Factor_Inverse);
		}
		BOOST_ASSERT(rs.front() == 0);
		BOOST_ASSERT(rs.back() == scoring_function::Cutoff);

		// Populate the scoring function task container.
		const size_t num_sf_tasks = ((XS_TYPE_SIZE + 1) * XS_TYPE_SIZE) >> 1;
		cnt.init(num_sf_tasks);
		for (size_t t0 =  0; t0 < XS_TYPE_SIZE; ++t0)
		for (size_t t1 = t0; t1 < XS_TYPE_SIZE; ++t1)
		{
			io.post([&,t0,t1]()
			{
				sf.precalculate(t0, t1, rs);
				cnt.increment();
			});
		}
		cnt.wait();
	}
	cout << '\n';

	cout << "Running " << num_mc_tasks << " Monte Carlo task" << ((num_mc_tasks == 1) ? "" : "s") << " per ligand\n";

	// Perform docking for each file in the ligand folder.
	log_engine log;
	cout.setf(ios::fixed, ios::floatfield);
	cout << "   Index        Ligand    pKd 1     2     3     4     5     6     7     8     9" << endl << setprecision(2);
	path input_ligand_path;
	size_t num_conformations; // Number of conformation to output.
	using namespace boost::filesystem;
	const directory_iterator end_dir_iter; // A default constructed directory_iterator acts as the end iterator.
	for (directory_iterator dir_iter(input_folder_path); dir_iter != end_dir_iter; ++dir_iter)
	{
		// Obtain a ligand.
		input_ligand_path = dir_iter->path();

		// Skip non-regular files such as folders.
		if (!is_regular_file(input_ligand_path)) continue;

		// Filter files with .pdbqt extension name.
		if (input_ligand_path.extension() != ".pdbqt") continue;

		const path output_ligand_path = output_folder_path / input_ligand_path.filename();

		// Parse the ligand.
		ligand lig(input_ligand_path);

		// Create grid maps on the fly if necessary.
		BOOST_ASSERT(atom_types_to_populate.empty());
		const vector<size_t> ligand_atom_types = lig.get_atom_types();
		const size_t num_ligand_atom_types = ligand_atom_types.size();
		for (size_t i = 0; i < num_ligand_atom_types; ++i)
		{
			const size_t t = ligand_atom_types[i];
			BOOST_ASSERT(t < XS_TYPE_SIZE);
			array3d<fl>& grid_map = grid_maps[t];
			if (grid_map.initialized()) continue; // The grid map of XScore atom type t has already been populated.
			grid_map.resize(b.num_probes); // An exception may be thrown in case memory is exhausted.
			atom_types_to_populate.push_back(t);  // The grid map of XScore atom type t has not been populated and should be populated now.
		}
		const size_t num_atom_types_to_populate = atom_types_to_populate.size();
		if (num_atom_types_to_populate)
		{
			// Populate the grid map task container.
			cnt.init(num_gm_tasks);
			for (size_t x = 0; x < num_gm_tasks; ++x)
			{
				io.post([&,x]()
				{
					grid_map_task(grid_maps, atom_types_to_populate, x, sf, b, rec);
					cnt.increment();
				});
			}

			// Block until all the grid map tasks are completed.
			cnt.wait();
			atom_types_to_populate.clear();
		}

		// Dump the ligand file stem.
		string stem = input_ligand_path.stem().string();
		cout << setw(8) << log.size() + 1 << setw(14) << stem << "   " << flush;

		// Populate the Monte Carlo task container.
		cnt.init(num_mc_tasks);
		for (size_t i = 0; i < num_mc_tasks; ++i)
		{
			BOOST_ASSERT(result_containers[i].empty());
			io.post([&,i]()
			{
				monte_carlo_task(result_containers[i], lig, eng(), sf, b, grid_maps);
				cnt.increment();
			});
		}
		cnt.wait();

		// Merge results from all the tasks into one single result container.
		BOOST_ASSERT(results.empty());
		const fl required_square_error = static_cast<fl>(4 * lig.num_heavy_atoms); // Ligands with RMSD < 2.0 will be clustered into the same cluster.
		for (size_t i = 0; i < num_mc_tasks; ++i)
		{
			ptr_vector<result>& task_results = result_containers[i];
			const size_t num_task_results = task_results.size();
			for (size_t j = 0; j < num_task_results; ++j)
			{
				add_to_result_container(results, static_cast<result&&>(task_results[j]), required_square_error);
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
		const fl best_result_intra_e = best_result.e - best_result.f;
		for (size_t i = 0; i < num_results; ++i)
		{
			results[i].e_nd = (results[i].e - best_result_intra_e) * lig.flexibility_penalty_factor;
		}

		// Write models to file.
		lig.write_models(output_ligand_path, results, num_results, b, grid_maps);

		// Display the free energies of the top 4 conformations.
		const size_t num_energies = min<size_t>(num_results, 9);
		for (size_t i = 0; i < num_energies; ++i)
		{
			cout << setw(6) << results[i].e_nd;
		}

		// Add a log record.
		vector<fl> affinities(num_results);
		for (size_t i = 0; i < num_results; ++i)
		{
			affinities[i] = results[i].e_nd;
		}
		log.push_back(new log_record(move(stem), move(affinities)));

		cout << endl;

		// Clear the results of the current ligand.
		results.clear();
	}

	// Sort and write ligand log records to the log file.
	if (log.empty()) return 0;
	cout << "Writing log records of " << log.size() << " ligands to " << log_path << endl;
	log.sort();
	log.write(log_path);
}
