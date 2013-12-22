#include <cstdlib>
#include <fstream>
#include "source.hpp"

source::source()
{
	ifstream ifs(getenv("idock_cl"));
	insert(begin(), istreambuf_iterator<char>(ifs), istreambuf_iterator<char>());
}
