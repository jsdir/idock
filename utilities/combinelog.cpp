#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

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
	string slice;
	size_t comma;
	double energy;

	explicit record(string&& line_, const string& slice) : line(static_cast<string&&>(line_)), slice(slice), comma(line.find_first_of(',', 12)), energy(lexical_cast<double>(line.substr(comma + 3, line.find_first_of(',', comma + 8) - comma - 3))) {}

	record(const record&) = default;
	record(record&&) = default;
	record& operator=(const record&) = default;
	record& operator=(record&&) = default;
};

inline bool operator<(const record& a, const record& b)
{
	return a.energy < b.energy;
}

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

class property
{
public:
	string line;
	string id;

	explicit property(string&& line_) : line(static_cast<string&&>(line_)), id(line.substr(4, 8)) {}
};

inline size_t binary(const vector<property>& properties, const string& id)
{
	size_t s = 0;
	size_t e = properties.size();
	while (s + 1 < e)
	{
		const size_t mid = (s + e) / 2;
		if (id < properties[mid].id)
		{
			e = mid;
		}
		else
		{
			s = mid;
		}
	}
	return s;
}

int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		std::cout << "combinelog slices_folder prefix 16_prop_350.xls out.csv\n";
		return 1;
	}

	const path slices = argv[1];
	const string prefix = argv[2];
	const path prop = argv[3];
	const path output_csv = argv[4];

	string line;
	line.reserve(600);

	std::cout << "Reading log.csv's." << std::endl;
	vector<record> records;
	records.reserve(80000);
	const directory_iterator end_dir_iter;
	for (directory_iterator dir_iter(slices); dir_iter != end_dir_iter; ++dir_iter)
	{
		if (!is_directory(dir_iter->status())) continue; // Find example directories.
		const path slice_path = dir_iter->path();
		const string fullslice = slice_path.filename().string();
		if (!starts_with(fullslice, prefix)) continue;
		std::cout << fullslice << '\n';
		const string slice = fullslice.substr(6);
		ifstream log(slice_path / "log.csv");
		getline(log, line); // Filter out header line.
		while (getline(log, line))
		{
			records.push_back(record(static_cast<string&&>(line), slice));
		}
	}

	const size_t num_records = records.size();
	std::cout << "Sorting " << num_records << " records." << std::endl;
	std::sort(records.begin(), records.end());

	std::cout << "Reading property xls file " << prop << '.' << std::endl;
	vector<property> properties;
	properties.reserve(7220928); // Number of molecules in 16_prop_350.xls
	ifstream xls(prop);
	getline(xls, line); // Filter out header line.
	while (getline(xls, line))
	{
		properties.push_back(property(static_cast<string&&>(line)));
	}
	xls.close();
	std::cout << properties.size() << " records in prop xls file." << std::endl;

	std::cout << "Writing combined csv to " << output_csv << std::endl;
	ofstream csv(output_csv);
	csv << "Slice,Ligand,Conf,FE1,HB1,FE2,HB2,FE3,HB3,FE4,HB4,FE5,HB5,FE6,HB6,FE7,HB7,FE8,HB8,FE9,HB9,HA,MWT,LogP,Desolv_apolar,Desolv_polar,HBD,HBA,tPSA,Charge,NRB,SMILES\n";
	size_t s, e, num_not_found = 0;
	for (size_t i = 0; i < num_records; ++i)
	{
		const record& r = records[i];
		const string id = r.line.substr(4, 8);
		csv << r.slice << ',' << id << ',' << r.line.substr(r.comma + 1);
		const property& p = properties[binary(properties, id)];
		if (id != p.id)
		{
			//std::cerr << id << ' '  << p.id;
			++num_not_found;
			csv << ",,,,,,,,,,,,,,,,,,,\n";
			continue;
		}
		s = 13;
		for (size_t j = 0; j < 10; ++j) // 9 properties, i.e. HA,MWT,LogP,Desolv_apolar,Desolv_polar,HBD,HBA,tPSA,Charge,NRB
		{
			e = p.line.find_first_of('\t', s + 1);
			csv << ',' << p.line.substr(s, e - s);
			s = e + 1;
		}
		csv << ',' << p.line.substr(s) << '\n'; // SMILES
	}
	csv.close();
	std::cout << num_not_found << " records in log.csv's but not in prop xls file." << std::endl;

	return 0;
}
