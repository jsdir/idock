#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "thread_pool.hpp"
#include "receptor.hpp"
#include "ligand.hpp"
#include "utility.hpp"
#include "kernel.hpp"
using namespace std;
using namespace boost::filesystem;

/// Represents a summary of docking results of a ligand.
class summary
{
public:
	const string stem;
	const float energy;
	explicit summary(const string& stem, const float energy) : stem(stem), energy(energy) {}
};

/// For sorting ptr_vector<summary>.
inline bool operator<(const summary& a, const summary& b)
{
	return a.energy < b.energy;
}

int main(int argc, char* argv[])
{
	path receptor_path, input_folder_path, output_folder_path, log_path;
	array<float, 3> center, size;
	size_t num_threads, seed, num_mc_tasks, num_generations, max_conformations;
	float granularity;

	// Parse program options in a try/catch block.
	try
	{
		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const path default_log_path = "log.csv";
		const size_t default_num_threads = thread::hardware_concurrency();
		const size_t default_seed = chrono::system_clock::now().time_since_epoch().count();
		const size_t default_num_mc_tasks = 256;
		const size_t default_num_generations = 100;
		const size_t default_max_conformations = 9;
		const float default_granularity = 0.15625f;

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
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "number of worker threads to use")
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("tasks", value<size_t>(&num_mc_tasks)->default_value(default_num_mc_tasks), "number of Monte Carlo tasks for global search")
			("generations", value<size_t>(&num_generations)->default_value(default_num_generations), "number of generations in BFGS")
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

	cout << "Creating a thread pool of " << num_threads << " worker thread" << (num_threads == 1 ? "" : "s") << endl;
	thread_pool tp(num_threads);

	cout << "Precalculating a scoring function of " << scoring_function::n << " atom types in parallel" << endl;
	scoring_function sf;
	for (size_t t2 = 0; t2 < sf.n; ++t2)
	for (size_t t1 = 0; t1 <=  t2; ++t1)
	{
		tp.push_back(packaged_task<int()>(bind(&scoring_function::precalculate, ref(sf), t1, t2)));
	}
	tp.sync();

	cout << "Parsing receptor " << receptor_path << endl;
	receptor rec(receptor_path, center, size, granularity);

	cout << "Setting up CUDA kernel with random seed " << seed << endl;
	kernel knl(sf.e.data(), sf.d.data(), sf.ns, sf.ne, rec.corner0.data(), rec.corner1.data(), rec.num_probes.data(), rec.granularity_inverse, num_mc_tasks, num_generations, seed);

	// Perform docking for each file in the ligand folder.
	boost::ptr_vector<summary> summaries;
	size_t num_ligands = 0; // Ligand counter.
	cout.setf(ios::fixed, ios::floatfield);
	cout << "Running " << num_mc_tasks << " Monte Carlo task" << (num_mc_tasks == 1 ? "" : "s") << " per ligand in parallel" << endl;
	cout << "   Index |          Ligand | kcal/mol" << endl << setprecision(2);
	const directory_iterator const_dir_iter; // A default constructed directory_iterator acts as the end iterator.
	for (directory_iterator dir_iter(input_folder_path); dir_iter != const_dir_iter; ++dir_iter)
	{
		// Parse the ligand.
		const path input_ligand_path = dir_iter->path();
		const ligand lig(input_ligand_path);

		// Output the ligand file stem.
		const string stem = input_ligand_path.stem().string();
		cout << setw(8) << ++num_ligands << " | " << setw(15) << stem << " | " << flush;

		// Create grid maps on the fly if necessary.
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
		if (xs.size())
		{
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
			for (size_t z = 0; z < rec.num_probes[2]; ++z)
			{
				tp.push_back(packaged_task<int()>(bind(&receptor::populate, ref(rec), cref(sf), cref(xs), z)));
			}
			tp.sync();
			knl.update(rec.maps, xs);
		}

		// Run the Monte Carlo tasks in parallel
		vector<float> ex;
		knl.launch(ex, lig.lig, lig.nv, lig.nf, lig.na, lig.np);

		// Cluster and save solutions to file.
		const path output_ligand_path = output_folder_path / input_ligand_path.filename();
		const float e = lig.save(output_ligand_path, ex, max_conformations, num_mc_tasks);
		cout << setw(8) << e << endl;
		summaries.push_back(new summary(stem, e));
	}

	// Sort and write ligand summary to the log file.
	cout << "Writing summary of " << num_ligands << " ligand" << (num_ligands == 1 ? "" : "s") << " to " << log_path << endl;
	summaries.sort();
	boost::filesystem::ofstream log(log_path);
	log.setf(ios::fixed, ios::floatfield);
	log << "Ligand,kcal/mol\n" << setprecision(2);
	for (size_t i = 0; i < num_ligands; ++i)
	{
		const summary& s = summaries[i];
		log << s.stem << ',' << s.energy << '\n';
	}
}
