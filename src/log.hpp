#pragma once
#ifndef IDOCK_LOG_HPP
#define IDOCK_LOG_HPP

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/filesystem/fstream.hpp>
using namespace std;
using namespace boost::filesystem;

/// Represents a log record of docking results of a ligand.
class log_record
{
public:
	const string stem;
	const vector<float> affinities;
	explicit log_record(string&& stem_, vector<float>&& affinities_) : stem(move(stem_)), affinities(move(affinities_)) {}
};

inline bool operator<(const log_record& r0, const log_record& r1)
{
	return r0.affinities.front() < r1.affinities.front();
}

class log_engine : public boost::ptr_vector<log_record>
{
public:
	using ptr_vector<log_record>::ptr_vector;
	/// Write ligand log records to the log file.
	void write(const path& log_path) const;
};

#endif
