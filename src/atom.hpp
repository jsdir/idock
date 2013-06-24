#pragma once
#ifndef IDOCK_ATOM_HPP
#define IDOCK_ATOM_HPP

#include <array>
#include <iomanip>
#include <boost/filesystem/fstream.hpp>
#include "utility.hpp"
using namespace std;
using namespace boost::filesystem;

class atom
{
private:
	static const size_t n = 29;
	static const array<string, n> ad_strings;
	static const array<float, n> ad_covalent_radii;
	static const array<size_t, n> ad_to_xs;
	static const array<size_t, n> ad_to_rf;
public:
	size_t serial;
	string name;
	array<float, 3> coord;
	size_t ad;
	size_t xs;
	size_t rf;

	explicit atom(const string& line) : serial(stoul(line.substr(6, 5))), name(line.substr(12, 4)), ad(find(ad_strings.cbegin(), ad_strings.cend(), line.substr(77, isspace(line[78]) ? 1 : 2)) - ad_strings.cbegin()), xs(ad_to_xs[ad]), rf(ad_to_rf[ad])
	{
		coord[0] = stof(line.substr(30, 8));
		coord[1] = stof(line.substr(38, 8));
		coord[2] = stof(line.substr(46, 8));
	}

	/// Returns the covalent radius of current AutoDock4 atom type.
	float covalent_radius() const
	{
		return ad_covalent_radii[ad];
	}

	/// Returns true if the AutoDock4 atom type is not supported.
	bool ad_unsupported() const
	{
		return ad == n;
	}

	/// Returns true if the XScore atom type is not supported.
	bool xs_unsupported() const
	{
		return xs == n;
	}

	/// Returns true if the RF-Score atom type is not supported.
	bool rf_unsupported() const
	{
		return rf == n;
	}

	/// Returns true if the atom is a nonpolar hydrogen atom.
	bool is_nonpolar_hydrogen() const
	{
		return ad == 0;
	}

	/// Returns true if the atom is a polar hydrogen atom.
	bool is_polar_hydrogen() const
	{
		return ad == 1;
	}

	/// Returns true if the atom is a hydrogen atom.
	bool is_hydrogen() const
	{
		return ad <= 1;
	}

	/// Returns true if the atom is a hetero atom, i.e. non-carbon heavy atom.
	bool is_hetero() const
	{
		return ad >= 4;
	}

	/// Returns true if the current atom is covalently bonded to a given atom.
	bool has_covalent_bond(const atom& a) const
	{
		return distance_sqr(coord, a.coord) < (covalent_radius() + a.covalent_radius()) * (covalent_radius() + a.covalent_radius());
	}

	/// For nitrogen and oxygen, revises the XScore atom type to make it a hydrogen bond donor.
	void donorize()
	{
		switch (xs)
		{
			case 2 : xs = 3; break; // Nitrogen, hydrogen bond donor.
			case 4 : xs = 5; break; // Nitrogen, both hydrogen bond donor and acceptor.
			case 6 : xs = 7; break; // Oxygen, both hydrogen bond donor and acceptor.
		}
	}

	/// For carbon, revises the XScore atom type to make it non-hydrophobic.
	void dehydrophobicize()
	{
		xs = 1; // Carbon, bonded to a hetero atom.
	}

	/// Output the atom in pdbqt format.
	void output(boost::filesystem::ofstream& ofs, const array<float, 3>& c) const
	{
		ofs << "ATOM  " << setw(5) << serial << ' ' << name << "              " << setw(8) << c[0] << setw(8) << c[1] << setw(8) << c[2] << "                " << "      " << ' ' << ad_strings[ad] << (ad_strings[ad].size() == 1 ? " " : "") << '\n';
	}
};

#endif
