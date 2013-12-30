#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>
using namespace std;

int main(int argc, char* argv[])
{
	const size_t d = 3;
	const array<size_t, d> p({ 30, 38, 46 });
	const array<char, d> c({ 'x', 'y', 'z' });
	vector<float> mn(d, numeric_limits<float>::max());
	vector<float> mx(d, numeric_limits<float>::lowest());
	for (string line; getline(cin, line);)
	{
		if (line[0] != 'A') continue;
		for (size_t i = 0; i < d; ++i)
		{
			const float v = stof(line.substr(p[i], 8));
			mn[i] = min<float>(mn[i], v);
			mx[i] = max<float>(mx[i], v);
		}
	}
	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(3);
	for (size_t i = 0; i < d; ++i)
	{
		cout << "center_" << c[i] << '=' << (mx[i] + mn[i]) * 0.5 << endl;
	}
	for (size_t i = 0; i < d; ++i)
	{
		cout << "size_"   << c[i] << '=' << (mx[i] - mn[i]) * 1.5 << endl;
	}
}