#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		cout << "filter id.txt < lines" << endl;
		return 1;
	}
	string id, line;
	for (ifstream ifs(argv[1]); getline(ifs, id);)
	{
		while (getline(cin, line))
		{
			if (line.substr(0, id.length()) == id)
			{
				cout << line << endl;
				break;
			}
		}
	}
}
