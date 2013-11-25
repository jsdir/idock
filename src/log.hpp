#pragma once
#ifndef IDOCK_LOG_HPP
#define IDOCK_LOG_HPP

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/filesystem/fstream.hpp>
using namespace boost;

/// Represents a log record of docking results of a ligand.
class log_record
{
public:
	const string stem;
	const vector<float> affinities;
	explicit log_record(string&& stem_, vector<float>&& affinities_) : stem(move(stem_)), affinities(move(affinities_)) {}
};

class log_engine : public ptr_vector<log_record>
{
public:
	using ptr_vector<log_record>::ptr_vector;
	/// Sort and write ligand log records to the log file.
	void write(const path& log_path);
};

#endif
