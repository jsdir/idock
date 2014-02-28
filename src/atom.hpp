#pragma once
#ifndef IDOCK_ATOM_HPP
#define IDOCK_ATOM_HPP

#include <array>
#include <string>
#include "vec3.hpp"
using std::string;

/// Represents an atom by very simple fields.
class atom
{
private:
	static const size_t n = 31; //!< Number of AutoDock4 atom types.
	static const std::array<string, n> ad_strings; //!< AutoDock4 atom type strings, e.g. H, HD, C, A.
	static const std::array<float, n> ad_covalent_radii; //!< Covalent radii of AutoDock4 atom types.
	static const std::array<size_t, n> ad_to_xs; //!< AutoDock4 to XScore atom type conversion.
public:
	size_t serial; ///< Serial number.
	string name; ///< Atom name;
	vec3 coordinate; ///< 3D coordinate.
	size_t ad; ///< AutoDock4 atom type.
	size_t xs; ///< XScore atom type.

	//! Constructs an atom from an ATOM/HETATM line in PDBQT format.
	explicit atom(const string& line);

	//! Returns true if the AutoDock4 atom type is not supported.
	bool ad_unsupported() const;

	/// Returns true if the atom is nonpolar hydrogen.
	bool is_nonpolar_hydrogen() const;

	/// Returns true if the atom is polar hydrogen.
	bool is_polar_hydrogen() const;

	/// Returns true if the atom is hydrogen.
	bool is_hydrogen() const;

	/// Returns true if the atom is a hetero atom, i.e. non-carbon heavy atom.
	bool is_hetero() const;

	/// Returns the covalent radius of current AutoDock4 atom type.
	double covalent_radius() const;

	/// Returns true if the current atom is covalently bonded to a given atom.
	bool is_neighbor(const atom& a) const;

	/// For nitrogen and oxygen, revises the XScore atom type to make it a hydrogen bond donor.
	void donorize();

	/// For carbon, revises the XScore atom type to make it non-hydrophobic.
	void dehydrophobicize();
};

#endif
