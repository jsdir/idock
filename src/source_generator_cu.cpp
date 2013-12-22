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
source<T>::source() : vector<T>({" << hex << setfill('0') << uppercase;
	for (char c; cin.get(c);)
	{
		cout << "0x" << setw(2) << static_cast<unsigned int>(static_cast<unsigned char>(c)) << ',';
	}
	cout << "})\n\
{\n\
}\n\
\n\
template class source<unsigned char>;\n\
";
}
