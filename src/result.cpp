#include "array.hpp"
#include "result.hpp"

//! Clusters a result into an existing result set with a minimum RMSD requirement.
// TODO: Consider using double linked list std::list<> to store results because of frequent insertions and deletions.
void result::push(ptr_vector<result>& results, result&& r, const double required_square_error)
{
	// If this is the first result, simply save it.
	if (results.empty())
	{
		results.push_back(new result(static_cast<result&&>(r)));
		return;
	}

	// If the container is not empty, find in a coordinate that is closest to the given newly found r.coordinate.
	size_t index = 0;
	double best_square_error = distance_sqr(r.heavy_atoms, results.front().heavy_atoms);
	for (size_t i = 1; i < results.size(); ++i)
	{
		const double this_square_error = distance_sqr(r.heavy_atoms, results[i].heavy_atoms);
		if (this_square_error < best_square_error)
		{
			index = i;
			best_square_error = this_square_error;
		}
	}

	if (best_square_error < required_square_error) // The result r is very close to results[index].
	{
		if (r.e < results[index].e) // r is better than results[index], so substitute r for results[index].
		{
			results.replace(index, new result(static_cast<result&&>(r)));
		}
	}
	else // Cannot find in results a result that is similar to r.
	{
		if (results.size() < results.capacity())
		{
			results.push_back(new result(static_cast<result&&>(r)));
		}
		else // Now the container is full.
		{
			if (r.e < results.back().e) // If r is better than the worst one, then replace it.
			{
				results.replace(results.size() - 1, new result(static_cast<result&&>(r)));
			}
		}
	}
	results.sort();
}
