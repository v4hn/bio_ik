#include "evolution.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

int seed = time(NULL);
boost::random::mt19937 rng(seed);

Evolution::Evolution(int size, int elites, int dimensionality, Dimension* &dimensions, FitnessFunction ff) {
	Size = size;
	Elites = elites;
	Dimensionality = dimensionality;
	Dimensions = dimensions;
	FF = ff;

	Population = new Individual[Size];
	Offspring = new Individual[Size];

	for(int i=0; i<Size; i++) {
		Population[i].Genes = new double[Dimensionality];
		Population[i].Gradient = new double[Dimensionality];
		Offspring[i].Genes = new double[Dimensionality];
		Offspring[i].Gradient = new double[Dimensionality];
		for(int j=0; j<Dimensionality; j++) {
			Population[i].Genes[j] = GetRandomValue(Dimensions[j].Min, Dimensions[j].Max);
			Population[i].Gradient[j] = 0.0;
			Offspring[i].Genes[j] = 0.0;
			Offspring[i].Gradient[j] = 0.0;
		}
		Population[i].Fitness = FF(Population[i].Genes, Dimensionality);
	}
	
	SortByFitness();
	ComputeExtinctions();

	Prototype.Genes = new double[Dimensionality];
	Prototype.Gradient = new double[Dimensionality];
	for(int i=0; i<Dimensionality; i++) {
		Prototype.Genes[i] = Population[0].Genes[i];
		Prototype.Gradient[i] = Population[0].Genes[i];
		Prototype.Extinction = Population[0].Extinction;
		Prototype.Fitness = Population[0].Fitness;
	}

	Initialized = true;
}

Evolution::Evolution(int size, int elites, int dimensionality, Dimension* &dimensions, FitnessFunction ff, vector<double> bias) {
	Size = size;
	Elites = elites;
	Dimensionality = dimensionality;
	Dimensions = dimensions;
	FF = ff;

	Population = new Individual[Size];
	Offspring = new Individual[Size];

	Population[0].Genes = new double[Dimensionality];
	Population[0].Gradient = new double[Dimensionality];
	for(int i=0; i<Dimensionality; i++) {
		Population[0].Genes[i] = bias[i];
		Population[0].Gradient[i] = 0.0;
	}
	Population[0].Fitness = FF(Population[0].Genes, Dimensionality);

	for(int i=1; i<Size; i++) {
		Population[i].Genes = new double[Dimensionality];
		Population[i].Gradient = new double[Dimensionality];
		for(int j=0; j<Dimensionality; j++) {
			Population[i].Genes[j] = GetRandomValue(Dimensions[j].Min, Dimensions[j].Max);
			Population[i].Gradient[j] = 0.0;
		}
		Population[i].Fitness = FF(Population[i].Genes, Dimensionality);
	}

	for(int i=0; i<Size; i++) {
		Offspring[i].Genes = new double[Dimensionality];
		Offspring[i].Gradient = new double[Dimensionality];
		for(int j=0; j<Dimensionality; j++) {
			Offspring[i].Genes[j] = 0.0;
			Offspring[i].Gradient[j] = 0.0;
		}
	}
	
	SortByFitness();
	ComputeExtinctions();

	Prototype.Genes = new double[Dimensionality];
	Prototype.Gradient = new double[Dimensionality];
	for(int i=0; i<Dimensionality; i++) {
		Prototype.Genes[i] = Population[0].Genes[i];
		Prototype.Gradient[i] = Population[0].Genes[i];
		Prototype.Extinction = Population[0].Extinction;
		Prototype.Fitness = Population[0].Fitness;
	}
}

Evolution::~Evolution() {
	for(int i=0; i<Size; i++) {
		delete[] Population[i].Genes;
		delete[] Population[i].Gradient;
	}
	delete[] Population;
	for(int i=0; i<Size; i++) {
		delete[] Offspring[i].Genes;
		delete[] Offspring[i].Gradient;
	}
	delete[] Offspring;

	delete[] Dimensions;
}

int Evolution::GetSize() {
	return Size;
}

int Evolution::GetDimensionality() {
	return Dimensionality;
}

void Evolution::Evolve() {
	//Let elites survive by performing a heuristic exploitation
	for(int i=0; i<Elites; i++) {
		Survive(i);
	}
	
	//Create mating pool
	vector<Individual*> pool;
	for(int i=0; i<Size; i++) {
		pool.push_back(&Population[i]);
	}
	
	//Evolve offspring
	for(int i=Elites; i<Size; i++) {
		if(pool.size() > 0) {
			Individual &parentA = Select(pool);
			Individual &parentB = Select(pool);

			Reproduce(i, parentA, parentB);

			if(Offspring[i].Fitness < parentA.Fitness) {
				//cout << "Removed A" << endl;
				pool.erase(remove(pool.begin(), pool.end(), &parentA), pool.end());
			}
			if(Offspring[i].Fitness < parentB.Fitness) {
				//cout << "Removed B" << endl;
				pool.erase(remove(pool.begin(), pool.end(), &parentB), pool.end());
			}
		} else {
			//cout << "Rerolling" << endl;
			Reroll(i);
		}
	}
	
	//Assign population
	swap(Population, Offspring);
	SortByFitness();

	//Assign extinction factors
	ComputeExtinctions();

	UpdatePrototype(Population[0]);

	//Check wipeout
	//if(CheckWipeout()) {
	//	cout << "WIPE SUGGESTED" <<endl;
	//}
}

void Evolution::UpdateFitnessFunction(FitnessFunction ff) {
	FF = ff;
}

Individual &Evolution::Select(vector<Individual*> &pool) {
	double* probabilities = new double[pool.size()];
	double rankSum = 0.0;
	for(int i=0; i<pool.size(); i++) {
		rankSum += (i+1);
	}
	for(int i=0; i<pool.size(); i++) {
		probabilities[i] = (pool.size()-i)/rankSum;
	}
	Individual &individual = *pool[GetRandomWeightedIndex(probabilities, pool.size())];
	delete[] probabilities;
	return individual;
}

void Evolution::Survive(int index) {
	Individual &survivor = Population[index];
	Individual &offspring = Offspring[index];

	for(int i=0; i<Dimensionality; i++) {
		offspring.Genes[i] = survivor.Genes[i];
		offspring.Gradient[i] = survivor.Gradient[i];
	}

	double fitnessSum = 0.0;
	for(int i=0; i<Dimensionality; i++) {
		double fitness = FF(offspring.Genes, Dimensionality);

		double inc = Clip(offspring.Genes[i] + GetRandomValue(0.0, fitness), Dimensions[i]) - offspring.Genes[i];
		offspring.Genes[i] += inc;
		double incFitness = FF(offspring.Genes, Dimensionality);
		offspring.Genes[i] -= inc;

		double dec = Clip(offspring.Genes[i] - GetRandomValue(0.0, fitness), Dimensions[i]) - offspring.Genes[i];
		offspring.Genes[i] += dec;
		double decFitness = FF(offspring.Genes, Dimensionality);
		offspring.Genes[i] -= dec;

		if(incFitness <= decFitness && incFitness < fitness) {
			offspring.Genes[i] += inc;
			offspring.Gradient[i] += inc;
			fitness = incFitness;
		}
		if(decFitness <= incFitness && decFitness < fitness) {
			offspring.Genes[i] += dec;
			offspring.Gradient[i] += dec;
			fitness = decFitness;
		}

		fitnessSum += fitness;
	}

	offspring.Fitness = fitnessSum/(double)Dimensionality;
}

void Evolution::Reproduce(int index, Individual &parentA, Individual &parentB) {
	Individual &offspring = Offspring[index];

	//Recombination
	for(int i=0; i<Dimensionality; i++) {
		float weight = GetRandomValue(0.0, 1.0);
		offspring.Genes[i] = weight*(parentA.Genes[i]) + (1.0-weight)*(parentB.Genes[i]) + GetRandomValue(-1.0, 1.0)*parentA.Gradient[i] + GetRandomValue(-1.0, 1.0)*parentB.Gradient[i];
	}

	double genesTmp[Dimensionality];
	for(int i=0; i<Dimensionality; i++) {
		genesTmp[i] = offspring.Genes[i];
	}

	//Mutation
	for(int i=0; i<Dimensionality; i++) {
		if(GetRandomValue(0.0, 1.0) < GetMutationProbability(parentA, parentB)) {
			offspring.Genes[i] += GetRandomValue(-1.0, 1.0) * GetMutationStrength(parentA, parentB, Dimensions[i]);
		}
	}

	//Adoption
	for(int i=0; i<Dimensionality; i++) {
		float weight = GetRandomValue(0.0, 1.0);
		offspring.Genes[i] += 
			weight * GetRandomValue(0.0, 1.0) * 0.5f * (parentA.Genes[i] - offspring.Genes[i] + parentB.Genes[i] - offspring.Genes[i])
			+ (1.0-weight) * GetRandomValue(0.0, 1.0) * (Population[0].Genes[i] - offspring.Genes[i]);
	}

	//Clip and Compute Evolutionary Gradient which is the change within Mutation and Adoption
	for(int i=0; i<Dimensionality; i++) {
		Clip(offspring);
		offspring.Gradient[i] = offspring.Genes[i] - genesTmp[i];
	}

	//Fitness
	offspring.Fitness = FF(offspring.Genes, Dimensionality);
}

void Evolution::Reroll(int index) {
	Individual &offspring = Offspring[index];

	for(int i=0; i<Dimensionality; i++) {
		offspring.Genes[i] = GetRandomValue(Dimensions[i].Min, Dimensions[i].Max);
		offspring.Gradient[i] = 0.0; 
	}

	offspring.Fitness = FF(offspring.Genes, Dimensionality);
}

Individual* &Evolution::GetPopulation() {
	return Population;
}

Individual &Evolution::GetPrototype() {
	return Prototype;
}

void Evolution::Clip(Individual &individual) {
	for(int i=0; i<Dimensionality; i++) {
		individual.Genes[i] = Clip(individual.Genes[i], Dimensions[i]);
	}
}

double Evolution::Clip(double value, const Dimension &dimension) {
	if(value < dimension.Min) {
		return dimension.Min;
	} else if(value > dimension.Max) {
		return dimension.Max;
	} else {
		return value;
	}
}

void Evolution::ComputeExtinctions() {
	double min = Population[0].Fitness;
	double max = Population[Size-1].Fitness;
	for(int i=0; i<Size; i++) {
		double grading = (double)i/((double)Size-1);
		Population[i].Extinction = (Population[i].Fitness + min*(grading-1.0)) / max;
	}
}

bool Evolution::CheckWipeout() {
	//In this, the balanced fitness function with equal weights must be used!
	if(FF(Prototype.Genes, Dimensionality) < FF(Population[0].Genes, Dimensionality)) {
		double* values = new double[Dimensionality];
		for(int i=0; i<Dimensionality; i++) {
			values[i] = Population[0].Genes[i];
		}
		for(int i=0; i<Dimensionality; i++) {
			double fitness = FF(values, Dimensionality);

			double inc = Clip(values[i] + GetRandomValue(0.0, fitness), Dimensions[i]) - values[i];
			values[i] += inc;
			double incFitness = FF(values, Dimensionality);
			values[i] -= inc;

			double dec = Clip(values[i] - GetRandomValue(0.0, fitness), Dimensions[i]) - values[i];
			values[i] += dec;
			double decFitness = FF(values, Dimensionality);
			values[i] -= dec;

			if(incFitness < fitness || decFitness < fitness) {
				delete[] values;
				return false;
			}
		}
		delete[] values;
		return true;
	}
	return false;
}

void Evolution::UpdatePrototype(Individual &candidate) {
	if(FF(candidate.Genes, Dimensionality) < FF(Prototype.Genes, Dimensionality)) {
		for(int i=0; i<Dimensionality; i++) {
			Prototype.Genes[i] = candidate.Genes[i];
			Prototype.Gradient[i] = candidate.Gradient[i];
			Prototype.Extinction = candidate.Extinction;
			Prototype.Fitness = candidate.Fitness;
		}
	}
}

void Evolution::SortByFitness() {
	sort(&Population[0], &Population[Size], SurvivalOfTheFittest());
}

double Evolution::GetMutationProbability(Individual &parentA, Individual &parentB) {
	double extinction = 0.5 * (parentA.Extinction + parentB.Extinction);
	double inverse = 1.0/Dimensionality;
	return extinction * (1.0-inverse) + inverse;
}

double Evolution::GetMutationStrength(Individual &parentA, Individual &parentB, Dimension &dimension) {
	double extinction = 0.5 * (parentA.Extinction + parentB.Extinction);
	double span = (dimension.Max - dimension.Min);
	return span * extinction;
}

double Evolution::GetRandomValue(double min, double max) {
	if(min == max) {
		return min;
	}
	boost::random::uniform_real_distribution<> gen(min, max);
	return gen(rng);
}

int Evolution::GetRandomWeightedIndex(double* &probabilities, int size) {
	double weightSum = 0.0;
	for(int i=0; i<size; i++) {
		weightSum += probabilities[i];
	}
	double rVal = GetRandomValue(0.0, 1.0)*weightSum;
	for(int i=0; i<size; i++) {
		rVal -= probabilities[i];
		if(rVal <= 0.0) {
			return i;
		}
	}
	return size-1;
}