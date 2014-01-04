#pragma once
#ifndef IDOCK_ATOM_HPP
#define IDOCK_ATOM_HPP

#include <array>
#include <boost/filesystem/fstream.hpp>
using namespace std;

//! Represents an atom of either receptor or ligand.
class atom
{
private:
	static const size_t n = 31; //!< Number of AutoDock4 atom types.
	static const array<string, n> ad_strings; //!< AutoDock4 atom type strings, e.g. H, HD, C, A.
	static const array<float, n> ad_covalent_radii; //!< Covalent radii of AutoDock4 atom types.
	static const array<size_t, n> ad_to_xs; //!< AutoDock4 to XScore atom type conversion.
	static const array<size_t, n> ad_to_rf; //!< AutoDock4 to RF-Score atom type conversion.
public:
	size_t serial; //!< Atom serial.
	string name; //!< Atom name, 4 characters wide.
	array<float, 3> coord; //!< Coordinate.
	size_t ad; //!< AutoDock4 atom type.
	size_t xs; //!< XScore atom type.
	size_t rf; //!< RF-Score atom type.
	vector<atom> hydrogens; //!< Hydrogens connected to the current atom.

	//! Constructs an atom from an ATOM/HETATM line in PDBQT format.
	explicit atom(const string& line);

	//! Returns true if the AutoDock4 atom type is not supported.
	bool ad_unsupported() const;

	//! Returns true if the XScore atom type is not supported.
	bool xs_unsupported() const;

	//! Returns true if the RF-Score atom type is not supported.
	bool rf_unsupported() const;

	//! Returns true if the atom is a nonpolar hydrogen atom, which is connected to a carbon atom.
	bool is_nonpolar_hydrogen() const;

	//! Returns true if the atom is a polar hydrogen atom, which is connected to a heavy atom other than carbon.
	bool is_polar_hydrogen() const;

	//! Returns true if the atom is a hydrogen atom, either polar or nonpolar.
	bool is_hydrogen() const;

	//! Returns true if the atom is a hetero atom, i.e. non-carbon heavy atom.
	bool is_hetero() const;

	//! Revises the XScore atom type of a nitrogen or an oxygen to its corresponding variant of hydrogen bond donor.
	void donorize();

	//! Revises the XScore atom type of a carbon to its corresponding non-hydrophobic variant.
	void dehydrophobicize();

	//! Returns the covalent radius of the current AutoDock4 atom type.
	float covalent_radius() const;

	//! Returns true if the current atom is covalently bonded to a given atom, i.e. their distance is within a certain threshold which depends on their covalent radii.
	bool has_covalent_bond(const atom& a) const;

	//! Outputs an ATOM line in PDBQT format.
	void output(boost::filesystem::ofstream& ofs, const array<float, 3>& coord) const;
};

#endif
