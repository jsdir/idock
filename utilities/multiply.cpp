#include <iostream>
#include <iomanip>
using namespace std;

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		cout << "multiply multiplier\n";
		return 1;
	}
	const double multiplier = stod(argv[1]);
	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(2);
	for (string line; getline(cin, line);)
	{
		cout << multiplier * stod(line) << endl;
	}
}