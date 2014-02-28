#pragma once
#ifndef IDOCK_MONTE_CARLO_TASK_HPP
#define IDOCK_MONTE_CARLO_TASK_HPP

#include <random>
#include "ligand.hpp"

const size_t num_alphas = 5; //!< Number of alpha values for determining step size in BFGS

//! Task for running Monte Carlo Simulated Annealing algorithm to find local minimums of the scoring function.
//! A Monte Carlo task uses a seed to initialize its own random number generator.
//! It starts from a random initial conformation,
//! repeats a specified number of iterations,
//! uses precalculated alpha values for line search during BFGS local search,
//! clusters free energies and heavy atom coordinate vectors of the best conformations into results,
//! and sorts the results in the ascending order of free energies.
void monte_carlo_task(ptr_vector<result>& results, const ligand& lig, const size_t seed, const scoring_function& sf, const box& b, const receptor& rec);

#endif
