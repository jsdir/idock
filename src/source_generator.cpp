#include <iostream>
#include <iomanip>
using namespace std;

int main()
{
	cout <<
"\
#include \"source.hpp\"\n\
\n\
source::source() : vector<char>({";
	for (char c; cin.get(c);)
	{
		cout << static_cast<int>(c) << ',';
	}
	cout << "})\n\
{\n\
}\n\
";
}
