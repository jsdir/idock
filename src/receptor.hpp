#pragma once
#ifndef IDOCK_RECEPTOR_HPP
#define IDOCK_RECEPTOR_HPP

#include <boost/filesystem/path.hpp>
#include "atom.hpp"
#include "scoring_function.hpp"
using namespace boost::filesystem;

//! Represents a receptor.
class receptor
{
public:
	vector<atom> atoms; //!< Heavy atoms.
	const array<float, 3> center; //!< Box center.
	const array<float, 3> size; //!< 3D sizes of box.
	const array<float, 3> corner0; //!< Box boundary corner with smallest values of all the 3 dimensions.
	const array<float, 3> corner1; //!< Box boundary corner with largest values of all the 3 dimensions.
	const float granularity; //!< 1D size of grids.
	const float granularity_inverse; //!< 1 / granularity.
	const array<int, 3> num_probes; //!< Number of probes.
	const size_t num_probes_product; //!< Product of num_probes[0,1,2].
	const size_t map_bytes; //!< Number of bytes in a map.
	vector<vector<size_t>> p_offset; //!< Auxiliary precalculated constants to accelerate grid map creation.
	vector<vector<float>> maps; //!< Grid maps.

	//! Constructs a receptor by parsing a receptor file in PDBQT format.
	explicit receptor(const path& p, const array<float, 3>& center, const array<float, 3>& size, const float granularity);

	//! Precalculates auxiliary constants to accelerate grid map creation.
	void precalculate(const scoring_function& sf, const vector<size_t>& xs);

	//! Populates grid maps for certain atom types along X and Y dimensions for a given Z dimension value.
	void populate(const vector<size_t>& xs, const size_t z, const scoring_function& sf);
};

#endif
