#ifndef EVOLUTION_H
#define EVOLUTION_H

#include <vector>
#include <iostream>
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/function.hpp>

using namespace std;

struct Dimension {
	double Min;
	double Max;
};

struct Individual {
    double* Genes;
	double* Gradient;
	double Extinction;
	double Fitness;
};

struct SurvivalOfTheFittest {
    inline bool operator() (const Individual &a, const Individual &b) {
        return (a.Fitness < b.Fitness);
    }
};

typedef boost::function<double(double* input, int dimensionality)> FitnessFunction;

class Evolution {
public:
	Evolution(int size, int elites, int dimensionality, Dimension* &dimensions, FitnessFunction ff);
	Evolution(int size, int elites, int dimensionality, Dimension* &dimensions, FitnessFunction ff, vector<double> bias);
	~Evolution();

	bool Initialized;

	Individual &GetPrototype();
	Individual* &GetPopulation();
	void Evolve();

	void UpdateFitnessFunction(FitnessFunction ff);

	int GetSize();
	int GetDimensionality();

private:
	Individual Prototype;

	int Size;
	int Elites;
	Dimension* Dimensions;
	int Dimensionality;
	FitnessFunction FF;

	Individual* Population;
	Individual* Offspring;

	Individual &Select(vector<Individual*> &pool);
	
	void Survive(int index);
	void Reproduce(int index, Individual &parentA, Individual &parentB);
	void Reroll(int index);
	
	void Clip(Individual &individual);
	double Clip(double value, const Dimension &dimension);
	void ComputeExtinctions();
	bool CheckWipeout();
	void UpdatePrototype(Individual &candidate);
	void SortByFitness();
	double GetMutationProbability(Individual &parentA, Individual &parentB);
	double GetMutationStrength(Individual &parentA, Individual &parentB, Dimension &dimension);

	double GetRandomValue(double min, double max);
	int GetRandomWeightedIndex(double* &probabilities, int size);
};
 
#endif