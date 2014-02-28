#pragma once
#ifndef IDOCK_BOX_HPP
#define IDOCK_BOX_HPP

#include <array>
using namespace std;

//! Represents a search space of cubic shape.
class box
{
public:
	static const double Default_Partition_Granularity; //!< Default size of partitions.
	static const double Default_Partition_Granularity_Inverse; //!< 1 / Default_Partition_Granularity.

	const array<double, 3> center; //!< Box center.
	array<double, 3> span; //!< 3D sizes of box.
	array<double, 3> corner1; //!< Box boundary corner with smallest values of all the 3 dimensions.
	array<double, 3> corner2; //!< Box boundary corner with largest values of all the 3 dimensions.
	const double grid_granularity; //!< 1D size of grids.
	const double grid_granularity_inverse; //!< 1 / grid_granularity.
	const array<double, 3> grid_size; //!< 3D sizes of grids.
	const array<double, 3> grid_size_inverse; //!< (1, 1, 1) / grid_size.
	array<size_t, 3> num_grids; //!< Number of grids.
	array<size_t, 3> num_probes; //!< Number of probes.
	array<size_t, 3> num_partitions; //!< Number of partitions.
	array<double, 3> partition_size; //!< 3D sizes of partitions.
	array<double, 3> partition_size_inverse; //!< (1, 1, 1) / partition_size.

	//! Constructs a search space of cubic shape.
	//! @param center Box center.
	//! @param span_ Intended 3D sizes of box. It will be expanded to the nearest multiple of grid_granularity.
	//! @param grid_granularity 1D size of grids.
	box(const array<double, 3>& center, const array<double, 3>& span_, const double grid_granularity);

	//! Returns true if a coordinate is within current half-open-half-close box, i.e. [corner1, corner2).
	bool within(const array<double, 3>& coordinate) const;

	//! Returns true if the distance between a coordinate and the surface of a box determined by boundary corner1 and corner2 is within cutoff.
	double project_distance_sqr(const array<double, 3>& corner1, const array<double, 3>& corner2, const array<double, 3>& coordinate) const;

	//! Returns true if the distance between a coordinate and the surface of current box is within cutoff.
	double project_distance_sqr(const array<double, 3>& coordinate) const;

	//! Returns the coordinate of boundary corner1 of the grid at the given 3D index.
	array<double, 3> grid_corner1(const array<size_t, 3>& index) const;

	//! Returns the coordinate of boundary corner1 of the partition at the given 3D index.
	array<double, 3> partition_corner1(const array<size_t, 3>& index) const;

	//! Returns the index of the half-open-half-close grid containing the given coordinate.
	array<size_t, 3> grid_index(const array<double, 3>& coordinate) const;

	//! Returns the index of the half-open-half-close partition containing the given coordinate.
	array<size_t, 3> partition_index(const array<double, 3>& coordinate) const;

	size_t grid_index(const array<size_t, 3>& a) const;
	size_t partition_index(const array<size_t, 3>& a) const;
};

#endif
