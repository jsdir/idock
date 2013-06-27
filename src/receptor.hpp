#pragma once
#ifndef IDOCK_RECEPTOR_HPP
#define IDOCK_RECEPTOR_HPP

#include <boost/filesystem/path.hpp>
#include "atom.hpp"
#include "array3d.hpp"
#include "scoring_function.hpp"
using namespace boost::filesystem;

/// Represents a receptor.
class receptor
{
public:
	static const float partition_granularity; ///< Default size of partitions.
	static const float partition_granularity_inv; ///< 1 / partition_granularity.

	vector<atom> atoms; ///< Receptor atoms.

	const array<float, 3> center; ///< Box center.
	array<float, 3> size; ///< 3D sizes of box.
	array<float, 3> corner0; ///< Box boundary corner with smallest values of all the 3 dimensions.
	array<float, 3> corner1; ///< Box boundary corner with largest values of all the 3 dimensions.
	const float grid_granularity; ///< 1D size of grids.
	const float grid_granularity_inverse; ///< 1 / grid_granularity.
	const array<float, 3> grid_size; ///< 3D sizes of grids.
	const array<float, 3> grid_size_inverse; ///< (1, 1, 1) / grid_size.
	array<size_t, 3> num_grids; ///< Number of grids.
	array<size_t, 3> num_probes; ///< Number of probes.
	array<size_t, 3> num_partitions; ///< Number of partitions.
	array<float, 3> partition_size; ///< 3D sizes of partitions.
	array<float, 3> partition_size_inverse; ///< (1, 1, 1) / partition_size.
	vector<array3d<float>> grid_maps;

	/// Constructs a receptor by parsing a receptor file in pdbqt format.
	/// @param center Box center.
	/// @param span_ Intended 3D sizes of box. It will be expanded to the nearest multiple of grid_granularity.
	/// @param grid_granularity 1D size of grids.
	explicit receptor(const path& p, const array<float, 3>& center, const array<float, 3>& size_, const float grid_granularity);

	/// Returns true if a coordinate is within current half-open-half-close box, i.e. [corner0, corner1).
	bool within(const array<float, 3>& coordinate) const;

	/// Returns true if the distance between a coordinate and the surface of a box determined by boundary corner0 and corner1 is within cutoff.
	float project_distance_sqr(const array<float, 3>& corner0, const array<float, 3>& corner1, const array<float, 3>& coordinate) const;

	/// Returns true if the distance between a coordinate and the surface of current box is within cutoff.
	float project_distance_sqr(const array<float, 3>& coordinate) const;

	/// Returns the coordinate of boundary corner0 of the grid at the given 3D index.
	array<float, 3> grid_corner0(const array<size_t, 3>& index) const;

	/// Returns the coordinate of boundary corner0 of the partition at the given 3D index.
	array<float, 3> partition_corner0(const array<size_t, 3>& index) const;

	/// Returns the index of the half-open-half-close grid containing the given coordinate.
	array<size_t, 3> grid_index(const array<float, 3>& coordinate) const;

	/// Returns the index of the half-open-half-close partition containing the given coordinate.
	array<size_t, 3> partition_index(const array<float, 3>& coordinate) const;

	/// Task for populating grid maps for certain atom types along X and Y dimensions for an Z dimension value.
	int populate(const scoring_function& sf, const vector<size_t>& xs, const size_t z);
};

#endif
