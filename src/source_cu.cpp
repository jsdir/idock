#include <cstdlib>
#include <fstream>
#include "source_cu.hpp"

source::source()
{
	ifstream ifs(getenv("idock_cu"), ios::binary);
	insert(begin(), (istreambuf_iterator<char>(ifs)), istreambuf_iterator<char>());
}
