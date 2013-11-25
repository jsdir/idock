#include "log.hpp"

/// For sorting ptr_vector<log_record>.
inline bool operator<(const log_record& r0, const log_record& r1)
{
	return r0.affinities.front() < r1.affinities.front();
}

void log_engine::write(const path& log_path)
{
	sort();
	boost::filesystem::ofstream log(log_path);
	log.setf(ios::fixed, ios::floatfield);
	log << "Ligand";
	for (size_t i = 1; i <= max_conformations; ++i)
	{
		log << ",pKd" << i;
	}
	log << '\n' << setprecision(2);
	for (const log_record& s : *this)
	{
		log << s.stem;
		for (const float a : s.affinities)
		{
			log << ',' << a;
		}
		for (size_t i = s.affinities.size(); i < max_conformations; ++i)
		{
			log << ',';
		}
		log << '\n';
	}
}
