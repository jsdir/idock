#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>

using namespace std;
using namespace boost;
using namespace boost::filesystem;

typedef double fl;

/// Returns true if a string starts with another string.
inline bool starts_with(const string& str, const string& start)
{
	const size_t start_size = start.size();
	if (str.size() < start_size) return false;
	for (size_t i = 0; i < start_size; ++i)
	{
		if (str[i] != start[i]) return false;
	}
	return true;
}

template<typename T>
T right_cast(const string& str, const size_t i, const size_t j)
{
	const size_t start = str.find_first_not_of(' ', i - 1);
	return lexical_cast<T>(str.substr(start, j - start));
}

class summary
{
public:
	const path filename;
	const vector<fl> energies;
	explicit summary(const path& filename, const vector<fl>& energies) : filename(filename), energies(energies) {}
};

/// For sorting ptr_vector<summary>.
inline bool operator<(const summary& a, const summary& b)
{
	return a.energies.front() < b.energies.front();
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		cout << "pdbqt2csv pdbqt_folder\n";
		return 1;
	}

	ptr_vector<summary> summaries;
	summaries.reserve(10000);
	vector<fl> energies;
	energies.reserve(9);
	string line;
	using boost::filesystem::directory_iterator;
	const directory_iterator end_dir_iter; // A default constructed directory_iterator acts as the end iterator.
	for (directory_iterator dir_iter(argv[1]); dir_iter != end_dir_iter; ++dir_iter)
	{
		// Skip non-regular files such as folders.
		if (!boost::filesystem::is_regular_file(dir_iter->status())) continue;
		const path p = dir_iter->path();
		boost::filesystem::fstream in(p);
		while (getline(in, line))
		{
			if (starts_with(line, "REMARK VINA RESULT:"))
			{
				energies.push_back(right_cast<fl>(line, 20, 29));
			}
			else if (starts_with(line, "REMARK       NORMALIZED FREE ENERGY PREDICTED BY IDOCK:"))
			{
				energies.push_back(right_cast<fl>(line, 56, 63));
			}
		}
		summaries.push_back(new summary(canonical(p), energies));
		energies.clear();
	}
	summaries.sort();
	cout << "ligand,no. of conformations";
	for (size_t i = 1; i <= 9; ++i)
	{
		cout << ",free energy in kcal/mol of conformation " << i;
	}
	cout.setf(std::ios::fixed, std::ios::floatfield);
	cout << '\n' << std::setprecision(3);
	const size_t num_summaries = summaries.size();
	for (size_t i = 0; i < num_summaries; ++i)
	{
		const summary& s = summaries[i];
		const size_t num_conformations = s.energies.size();
		cout << s.filename << ',' << num_conformations;
		for (size_t j = 0; j < num_conformations; ++j)
		{
			cout << ',' << s.energies[j];
		}
		cout << '\n';
	}
	
	return 0;
}
