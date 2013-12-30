#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>

using std::string;
using boost::lexical_cast;
using boost::filesystem::path;
using boost::filesystem::ifstream;
using boost::filesystem::ofstream;

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

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "extractmodel models.pdbqt model\n";
		return 1;
	}

	const path models = argv[1];
	const size_t model = lexical_cast<size_t>(argv[2]);
	const string model_start = "MODEL        " + lexical_cast<string>(model);
	const string model_end = "ENDMDL";

	string line;
	ofstream out(models.filename());
	ifstream in(models);
	while (getline(in, line) && !starts_with(line, model_start));
	while (getline(in, line) && !starts_with(line, model_end))
	{
		out << line << '\n';
	}
	return 0;
}
