#include <iostream>
#include <iomanip>
using namespace std;

int main()
{
	cout <<
"\
#include \"source.hpp\"\n\
\n\
template <typename T>\n\
source<T>::source() : vector<T>({";
	for (char c; cin.get(c);)
	{
		cout << static_cast<unsigned int>(c) << ',';
	}
	cout << "})\n\
{\n\
}\n\
\n\
template class source<char>;\n\
";
}
