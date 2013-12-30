#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

using std::cout;
using std::string;
using std::vector;
using boost::lexical_cast;
using boost::filesystem::path;
using boost::filesystem::directory_iterator;
using boost::filesystem::is_directory;
using boost::filesystem::ifstream;
using boost::filesystem::ofstream;

class record
{
public:
	string line;
	double energy;
	bool hb;

	explicit record(string&& line_, const bool hb) : line(static_cast<string&&>(line_)), hb(hb)
	{
		const auto comma = line.find(',', 12);
		energy = lexical_cast<double>(line.substr(comma + 1, line.find(',', comma + 7) - comma - 1));
	}

	record(const record&) = default;
	record(record&&) = default;
	record& operator=(const record&) = default;
	record& operator=(record&&) = default;
};

inline bool operator<(const record& a, const record& b)
{
	return a.energy < b.energy;
}

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "combinelog logs_folder out.csv\n";
		return 1;
	}

	string line;
	line.reserve(600);

	std::cout << "Reading log.csv's." << std::endl;
	vector<record> records;
	records.reserve(7220836);
	const directory_iterator end_dir_iter;
	for (directory_iterator dir_iter(argv[1]); dir_iter != end_dir_iter; ++dir_iter)
	{
		const path log_path = dir_iter->path();
		std::cout << log_path << '\n';
		ifstream log(log_path);
		getline(log, line); // Header line.
		const auto hb = line.size() == 156;
		while (getline(log, line))
		{
			records.push_back(record(static_cast<string&&>(line), hb));
		}
	}

	const size_t num_records = records.size();
	std::cout << "Sorting " << num_records << " records." << std::endl;
	std::sort(records.begin(), records.end());

	std::cout << "Writing combined csv\n";
	ofstream csv(argv[2]);
	csv << "Slice,Ligand,Conf,FE1,HB1,FE2,HB2,FE3,HB3,FE4,HB4,FE5,HB5,FE6,HB6,FE7,HB7,FE8,HB8,FE9,HB9,MWT,LogP,Desolv_apolar,Desolv_polar,HBD,HBA,tPSA,Charge,NRB,SMILES\n";
	for (size_t i = 0; i < num_records; ++i)
	{
		const record& r = records[i];
		if (r.hb)
		{
			csv << r.line;
		}
		else
		{
			size_t s = r.line[12] == ',' ? 13 : 14;
			csv << r.line.substr(0, s);
			size_t e;
			for (size_t j = 0; j < 9; ++j) // 9 properties, i.e. MWT,LogP,Desolv_apolar,Desolv_polar,HBD,HBA,tPSA,Charge,NRB
			{
				e = r.line.find(',', s);
				csv << r.line.substr(s, e - s) << ",,";
				s = e + 1;
			}
			csv << r.line.substr(s);
		}
		csv << '\n';
	}
	csv.close();

	return 0;
}
