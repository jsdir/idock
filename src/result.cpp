#include "array.hpp"
#include "result.hpp"

//! Clusters a result into a result container with a minimum RMSD requirement.
void result::push(ptr_vector<result>& results, result&& r, const double required_square_error)
{
	// If this is the first result, simply save it.
	if (results.empty())
	{
		results.push_back(new result(move(r)));
		return;
	}

	// If the container is not empty, find a result to which r is the closest.
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

	// Now r is the closest to results[index]. Check if they are in the same cluster.
	if (best_square_error < required_square_error)
	{
		// They are in the same cluster and r is better than results[index], so substitute r for results[index].
		if (r.e < results[index].e)
		{
			results.replace(index, new result(move(r)));
		}
	}
	else // They are not in the same cluster, i.e. r itself forms a new cluster.
	{
		// Save this new cluster if the result container is not full yet.
		if (results.size() < results.capacity())
		{
			results.push_back(new result(move(r)));
		}
		else // Now the container is full.
		{
			// If r is better than the worst result, then substitute r for it.
			if (r.e < results.back().e)
			{
				results.replace(results.size() - 1, new result(move(r)));
			}
		}
	}
	results.sort();
}
