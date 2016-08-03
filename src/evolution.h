#ifndef EVOLUTION_H
#define EVOLUTION_H

#include <vector>
#include <iostream>
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <tf_conversions/tf_kdl.h>
#include <kdl/chain.hpp>

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

class Evolution {
public:
	Evolution(int populationSize, int elites, int dimensionality, Dimension* dimensions, const vector<double> &seed, const geometry_msgs::Pose &target, const KDL::Chain &chain, const double &chainLength);
	~Evolution();

	double* &GetSolution();
	Individual* &GetPopulation();
	void Evolve();
	void Terminate();

	int GetPopulationSize();
	int GetDimensionality();
	double GetSolutionFitness();

private:
	geometry_msgs::Pose Target;
	KDL::Chain Chain;
	double ChainLength;
	int JointCount;
	int SegmentCount;
	
	geometry_msgs::Pose BasePose, EEPose;

	vector<double> Seed;

	int PopulationSize;
	int Elites;
	const Dimension* Dimensions;
	int Dimensionality;

	double* Solution;
	double SolutionFitness;

	Individual* Population;
	Individual* Offspring;

	void Initialize();

	Individual &Select(vector<Individual*> &pool);
	
	void Survive(int &index);
	void Reproduce(int &index, Individual &parentA, Individual &parentB);
	void Reroll(int &index);

	double ComputeFitness(double* &genes);
	double ComputeFitness(KDL::Frame &frame);
	double ComputeBalancedFitness(double* &genes);
	double ComputeBalancedFitness(KDL::Frame &frame);
	double GetPoseFitness(double balance);
	void ComputeFK(double* &values);

	//void Clip(Individual &individual);
	double Clip(double value, const Dimension &dimension);
	void ComputeExtinctions();
	bool CheckWipeout();
	bool UpdateSolution(Individual &candidate);
	void SortByFitness();
	double GetMutationProbability(Individual &parentA, Individual &parentB);
	double GetMutationStrength(Individual &parentA, Individual &parentB, const Dimension &dimension);

	double GetRandomValue(double min, double max);
	int GetRandomWeightedIndex(double* &probabilities, int size);
	double GetAngleDifference(double& q1x, double& q1y, double& q1z, double& q1w, double& q2x, double& q2y, double& q2z, double& q2w);
};
 
#endif