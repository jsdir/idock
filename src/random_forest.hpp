#pragma once
#ifndef IDOCK_RANDOM_FOREST_HPP
#define IDOCK_RANDOM_FOREST_HPP

#include <vector>
#include <array>
#include <random>
#include <mutex>
using namespace std;

//! Represents a node in a tree.
class node
{
public:
	vector<size_t> samples; //!< Node samples.
	float y; //!< Average of y values of node samples.
	float p; //!< Node purity, measured as either y * y * nSamples or sum * sum / nSamples.
	size_t var; //!< Variable used for node split.
	float val; //!< Value used for node split.
	array<size_t, 2> children; //!< Two child nodes.
};

//! Represents a tree in a forest.
class tree : public vector<node>
{
public:
	static const size_t nv = 41; //!< Number of variables.

	//! Trains an empty tree from bootstrap samples.
	void train(const size_t mtry, const function<float()> u01);

	//! Predicts the y value of the given sample x.
	float operator()(const array<float, nv>& x) const;

	//! Clears node samples to save memory.
	void clear();
private:
	static const size_t ns = 2959; //!< Number of training samples.
	static const array<array<float, nv>, ns> x; //!< Training samples.
	static const array<float, ns> y; //!< y values of training samples.
};

//! Represents a random forest.
class forest : public vector<tree>
{
public:
	//! Constructs a random forest of a number of empty trees.
	forest(const size_t nt, const size_t seed);

	//! Predicts the y value of the given sample x.
	float operator()(const array<float, tree::nv>& x) const;

	//! Clears node samples to save memory.
	void clear();

	//! Returns a random value from uniform distribution in [0, 1] in a thread safe manner.
	const function<float()> u01_s;
private:
	mt19937_64 rng;
	uniform_real_distribution<float> uniform_01;
	mutable mutex m;
};

#endif
