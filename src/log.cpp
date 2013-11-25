#include <iomanip>
#include "log.hpp"

void log_engine::write(const path& log_path) const
{
	const size_t max_conformations = front().affinities.capacity();
	boost::filesystem::ofstream log(log_path);
	log.setf(ios::fixed, ios::floatfield);
	log << "Ligand";
	for (size_t i = 1; i <= max_conformations; ++i)
	{
		log << ",pKd" << i;
	}
	log << '\n' << setprecision(2);
	for (const auto& r : *this)
	{
		log << r.stem;
		for (const float a : r.affinities)
		{
			log << ',' << a;
		}
		for (size_t i = r.affinities.size(); i < max_conformations; ++i)
		{
			log << ',';
		}
		log << '\n';
	}
}
