#include <cstdlib>
#include <fstream>
#include "source.hpp"

template <typename T>
source<T>::source()
{
	ifstream ifs(getenv("idock_cu"), ios::binary);
	this->insert(this->begin(), istreambuf_iterator<char>(ifs), istreambuf_iterator<char>());
}

template class source<unsigned char>;
