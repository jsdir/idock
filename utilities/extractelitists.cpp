#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>

using std::string;
using boost::lexical_cast;
using boost::filesystem::path;
using boost::filesystem::ifstream;

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "extractelitists examples_folder num_elitists\n";
		return 1;
	}

	const path examples = argv[1];
	const size_t num_elitists = lexical_cast<size_t>(argv[2]);

	const string output = "output";
	const string pdbqtext = ".pdbqt";
	string line;
	line.reserve(600);

	ifstream log(examples / "log.csv");
	getline(log, line); // Filter out header line.
	for (size_t i = 0; i < num_elitists; ++i)
	{
		getline(log, line);
		const size_t comma = line.find_first_of(',', 1);
		const path ligand = "ZINC" + line.substr(comma + 1, 8) + pdbqtext;
		if (exists(ligand)) continue; // If this ligand has been extracted, no action is needed.
		const string slice = "16_p0." + line.substr(0, comma);
		copy(examples / slice / output / ligand, ligand); // Copy the elite ligand to the current working directory.
	}
	return 0;
}
