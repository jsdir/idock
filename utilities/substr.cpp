#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		cout << "substr pos len\n";
		return 1;
	}
	const size_t pos = stoul(argv[1]);
	const size_t len = stoul(argv[2]);
	const size_t end = pos + len;
	for (string line; getline(cin, line);)
	{
		if (line.size() < end) continue;
		cout << line.substr(pos, len) << endl;
	}
}