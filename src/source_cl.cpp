#include <cstdlib>
#include <fstream>
#include "source_cl.hpp"

source::source()
{
	ifstream ifs(getenv("idock_cl"));
	insert(begin(), (istreambuf_iterator<char>(ifs)), istreambuf_iterator<char>());
}
