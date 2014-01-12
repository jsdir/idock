#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

int main(int argc, char* argv[])
{
	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(8);
	string p, q;
	for (ifstream ifsp(argv[1]), ifsq(argv[2]); getline(ifsp, p) && getline(ifsq, q);)
	{
		cout << stod(p) / stod(q) << endl;
	}
}
