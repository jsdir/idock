#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
using namespace std;

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		cout << "rmsd reference.pdbqt docked.pdbqt\n";
		return 1;
	}

	string line;
	vector<double> ref_x, ref_y, ref_z;
	for (ifstream ifs(argv[1]); getline(ifs, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "ATOM  " || record == "HETATM")
		{
			const string element = line.substr(77, 2);
			if (!(element == "H " || element == "HD"))
			{
				ref_x.push_back(stod(line.substr(30, 8)));
				ref_y.push_back(stod(line.substr(38, 8)));
				ref_z.push_back(stod(line.substr(46, 8)));
			}
		}
	}
	const double n_inv = 1.0 / ref_x.size();

	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(2);
	double se = 0;
	size_t i = 0;
	for (ifstream ifs(argv[2]); getline(ifs, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "ATOM  " || record == "HETATM")
		{
			const string element = line.substr(77, 2);
			if (!(element == "H " || element == "HD"))
			{
				const double dx = ref_x[i] - stod(line.substr(30, 8));
				const double dy = ref_y[i] - stod(line.substr(38, 8));
				const double dz = ref_z[i] - stod(line.substr(46, 8));
				se += dx * dx + dy * dy + dz * dz;
				++i;
			}
		}
		else if (record == "TORSDO")
		{
			cout << sqrt(se * n_inv) << endl;
			se = i = 0;
		}
	}
}
