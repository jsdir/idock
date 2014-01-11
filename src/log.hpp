#pragma once
#ifndef IDOCK_LOG_HPP
#define IDOCK_LOG_HPP

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/filesystem/path.hpp>
using namespace std;
using namespace boost::filesystem;

//! Represents a log record of docking result of a ligand.
class log_record
{
public:
	const string stem; //!< Stem of the ligand filename.
	const vector<double> affinities; //!< Predicted binding affinities of the ligand.

	//! Constructs a log record by moving the file stem and predicted binding affinities of a ligand.
	explicit log_record(string&& stem_, vector<double>&& affinities_) : stem(move(stem_)), affinities(move(affinities_)) {}
};

//! Compares two log records by their first predicted binding affinity.
inline bool operator<(const log_record& r0, const log_record& r1)
{
	return r0.affinities.front() < r1.affinities.front();
}

//! Represents a vector of log records.
class log_engine : public boost::ptr_vector<log_record>
{
public:
	//! Write ligand log records to the log file.
	void write(const path& log_path) const;
};

#endif
