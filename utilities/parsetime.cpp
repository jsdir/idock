#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>

using namespace std;
using namespace boost;
using namespace boost::filesystem;

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		cout << "parsetime examples_folder\n";
		return 1;
	}

	cout << "Receptor,Example,Program,Command being timed,User time (seconds),System time (seconds),Percent of CPU this job got,Elapsed (wall clock) time (h:mm:ss or m:ss),Average shared text size (kbytes),Average unshared data size (kbytes),Average stack size (kbytes),Average total siz (kbytes),Maximum resident set size (kbytes),Average resident set size (kbytes),Major (requiring I/O) page faults,Minor (reclaiming a frame) page faults,Voluntary context switches,Involuntary context switches,Swaps,File system inputs,File system outputs,Socket messages sent,Socket messages received,Signals delivered,Page size (bytes),Exit status\n";
	const string programs[2] = { "vina", "idock" };
	string line;
	using boost::filesystem::directory_iterator;
	const directory_iterator end_dir_iter;
	for (directory_iterator receptor_dir_iter(argv[1]); receptor_dir_iter != end_dir_iter; ++receptor_dir_iter)
	{
		// Find receptor directories.
		if (!boost::filesystem::is_directory(receptor_dir_iter->status())) continue;
		const path receptor = receptor_dir_iter->path();
		for (directory_iterator example_dir_iter(receptor); example_dir_iter != end_dir_iter; ++example_dir_iter)
		{
			// Find example directories.
			if (!boost::filesystem::is_directory(example_dir_iter->status())) continue;
			const path example = example_dir_iter->path();
			const string example_string = example.filename().string();
			if (example_string[0] != '1') continue; // Only consider 16_*
			for (size_t i = 0 ; i < 2; ++i)
			{
				const string& program = programs[i];
				cout << receptor.stem().string() << ',' << example_string << ',' << program;
				boost::filesystem::fstream in(example / (program + ".time"));
				getline(in, line);
				cout << ',' << line.substr(line.find_first_of(":") + 2); // e.g. Command being timed: "idock --config idock.cfg"
				while (getline(in, line))
				{
					cout << ',' << line.substr(line.find_last_of(" ") + 1);
				}
				cout << '\n';
			}
		}
	}
	return 0;
}
