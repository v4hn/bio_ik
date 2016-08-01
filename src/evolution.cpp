#include "evolution.h"

boost::random::mt19937 RNG(time(NULL));
boost::random::uniform_real_distribution<> ZeroOneUniform(0.0, 1.0);

Evolution::Evolution(int populationSize, int elites, int dimensionality, Dimension* dimensions, const vector<double> &seed, const geometry_msgs::Pose &target, const KDL::Chain &chain, const double &chainLength, const int &jointCount, const int &segmentCount) {
	Target = target;
	Chain = chain;
	ChainLength = chainLength;
	JointCount = jointCount;
	SegmentCount = segmentCount;

	PopulationSize = populationSize;
	Elites = elites;
	Dimensionality = dimensionality;
	Dimensions = dimensions;

	Population = new Individual[PopulationSize];
	Population[0].Genes = new double[Dimensionality];
	Population[0].Gradient = new double[Dimensionality];
	for(int i=0; i<Dimensionality; i++) {
		Population[0].Genes[i] = seed[i];
		Population[0].Gradient[i] = 0.0;
	}
	Population[0].Fitness = ComputeFitness(Population[0].Genes);
	for(int i=1; i<PopulationSize; i++) {
		Population[i].Genes = new double[Dimensionality];
		Population[i].Gradient = new double[Dimensionality];
		for(int j=0; j<Dimensionality; j++) {
			Population[i].Genes[j] = GetRandomValue(Dimensions[j].Min, Dimensions[j].Max);
			Population[i].Gradient[j] = 0.0;
		}
		Population[i].Fitness = ComputeFitness(Population[i].Genes);
	}

	Offspring = new Individual[PopulationSize];
	for(int i=0; i<PopulationSize; i++) {
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

	tf::poseKDLToMsg(Chain.getSegment(0).pose(0.0), BasePose);

	EvolutionFitness = ComputeBalancedFitness(Prototype.Genes);
}

Evolution::~Evolution() {
	for(int i=0; i<PopulationSize; i++) {
		delete[] Population[i].Genes;
		delete[] Population[i].Gradient;
	}
	delete[] Population;
	for(int i=0; i<PopulationSize; i++) {
		delete[] Offspring[i].Genes;
		delete[] Offspring[i].Gradient;
	}
	delete[] Offspring;

	//delete[] Dimensions;
}

int Evolution::GetPopulationSize() {
	return PopulationSize;
}

int Evolution::GetDimensionality() {
	return Dimensionality;
}

double Evolution::GetEvolutionFitness() {
	return EvolutionFitness;
}

void Evolution::Evolve() {
	//Let elites survive by performing a heuristic exploitation
	for(int i=0; i<Elites; i++) {
		Survive(i);
	}

	//Create mating pool
	vector<Individual*> pool;
	for(int i=0; i<PopulationSize; i++) {
		pool.push_back(&Population[i]);
	}
	
	//Evolve offspring
	for(int i=Elites; i<PopulationSize; i++) {
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

Individual &Evolution::Select(vector<Individual*> &pool) {
	int size = pool.size();
	double* probabilities = new double[size];
	double rankSum = 0.0;
	for(int i=0; i<size; i++) {
		rankSum += (i+1);
	}
	for(int i=0; i<size; i++) {
		probabilities[i] = (size-i)/rankSum;
	}
	Individual &individual = *pool[GetRandomWeightedIndex(probabilities, size)];
	delete[] probabilities;
	return individual;
}

void Evolution::Survive(int &index) {
	Individual &survivor = Population[index];
	Individual &offspring = Offspring[index];

	for(int i=0; i<Dimensionality; i++) {
		offspring.Genes[i] = survivor.Genes[i];
		offspring.Gradient[i] = survivor.Gradient[i];
	}

	//Perform Exploitation (Slow)
	/*
	double fitnessSum = 0.0;
	for(int i=0; i<Dimensionality; i++) {
		double fitness = ComputeFitness(offspring.Genes);

		double inc = Clip(offspring.Genes[i] + GetRandomValue(0.0, fitness), Dimensions[i]) - offspring.Genes[i];
		offspring.Genes[i] += inc;
		double incFitness = ComputeFitness(offspring.Genes);
		offspring.Genes[i] -= inc;

		double dec = Clip(offspring.Genes[i] - GetRandomValue(0.0, fitness), Dimensions[i]) - offspring.Genes[i];
		offspring.Genes[i] += dec;
		double decFitness = ComputeFitness(offspring.Genes);
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
	*/
	
	//Perform Exploitation (Fast)
	double fitnessSum = 0.0;
	
	KDL::Frame current, rest, result;
	KDL::Frame defFrame, incFrame, decFrame;
	int dimension;

	dimension = 0;
	for(int i=0; i<SegmentCount; i++) {
		if(Chain.getSegment(i).getJoint().getType() != KDL::Joint::None) {
			rest = rest*Chain.getSegment(i).pose(offspring.Genes[dimension]);
			dimension += 1;
		} else {
			rest = rest*Chain.getSegment(i).pose(0.0);
		}
	}
	
	dimension = 0;
	for(int i=0; i<SegmentCount; i++) {
		if(Chain.getSegment(i).getJoint().getType() != KDL::Joint::None) {
			result = current*rest;
			double fitness = ComputeFitness(result);

			defFrame = Chain.getSegment(i).pose(offspring.Genes[dimension]);
			rest = defFrame.Inverse()*rest;

			double inc = Clip(offspring.Genes[dimension] + GetRandomValue(0.0, fitness), Dimensions[dimension]) - offspring.Genes[dimension];
			incFrame = Chain.getSegment(i).pose(offspring.Genes[dimension]+inc);
			result = current*incFrame*rest;
			double incFitness = ComputeFitness(result);

			double dec = Clip(offspring.Genes[dimension] - GetRandomValue(0.0, fitness), Dimensions[dimension]) - offspring.Genes[dimension];
			decFrame = Chain.getSegment(i).pose(offspring.Genes[dimension]+dec);
			result = current*decFrame*rest;
			double decFitness = ComputeFitness(result);

			if(incFitness >= fitness && decFitness >= fitness) {
				current = current*defFrame;
			} else {
				if(incFitness < decFitness) {
					offspring.Genes[dimension] += inc;
					offspring.Gradient[dimension] += inc;
					fitness = incFitness;
					current = current*incFrame;
				} else {
					offspring.Genes[dimension] += dec;
					offspring.Gradient[dimension] += dec;
					fitness = decFitness;
					current = current*decFrame;
				}
			}

			fitnessSum += fitness;

			dimension += 1;
		} else {
			result = Chain.getSegment(i).pose(0.0);
			rest = result.Inverse()*rest;
			current = current*result;
		}
		offspring.Fitness = fitnessSum/(double)Dimensionality;
	}
}

void Evolution::Reproduce(int &index, Individual &parentA, Individual &parentB) {
	Individual &offspring = Offspring[index];

	//Recombination
	for(int i=0; i<Dimensionality; i++) {
		float weight = ZeroOneUniform(RNG);
		offspring.Genes[i] = weight*parentA.Genes[i] + (1.0-weight)*parentB.Genes[i] + GetRandomValue(-1.0, 1.0)*parentA.Gradient[i] + GetRandomValue(-1.0, 1.0)*parentB.Gradient[i];
	}

	double genesTmp[Dimensionality];
	for(int i=0; i<Dimensionality; i++) {
		genesTmp[i] = offspring.Genes[i];
	}

	//Mutation
	for(int i=0; i<Dimensionality; i++) {
		if(ZeroOneUniform(RNG) < GetMutationProbability(parentA, parentB)) {
			offspring.Genes[i] += GetRandomValue(-1.0, 1.0) * GetMutationStrength(parentA, parentB, Dimensions[i]);
		}
	}

	//Adoption
	for(int i=0; i<Dimensionality; i++) {
		float weight = ZeroOneUniform(RNG);
		offspring.Genes[i] += 
			weight * ZeroOneUniform(RNG) * 0.5f * (parentA.Genes[i] - offspring.Genes[i] + parentB.Genes[i] - offspring.Genes[i])
			+ (1.0-weight) * ZeroOneUniform(RNG) * (Population[0].Genes[i] - offspring.Genes[i]);
	}

	//Clip and Compute Evolutionary Gradient which is the change within Mutation and Adoption
	for(int i=0; i<Dimensionality; i++) {
		Clip(offspring);
		offspring.Gradient[i] = offspring.Genes[i] - genesTmp[i];
	}

	//Fitness
	offspring.Fitness = ComputeFitness(offspring.Genes);
}

void Evolution::Reroll(int &index) {
	Individual &offspring = Offspring[index];

	for(int i=0; i<Dimensionality; i++) {
		offspring.Genes[i] = GetRandomValue(Dimensions[i].Min, Dimensions[i].Max);
		offspring.Gradient[i] = 0.0; 
	}

	offspring.Fitness = ComputeFitness(offspring.Genes);
}

double Evolution::ComputeFitness(double* &genes) {
	ComputeFK(genes);
	
    //Position Error
    double diffX = EEPose.position.x - Target.position.x;
    double diffY = EEPose.position.y - Target.position.y;
    double diffZ = EEPose.position.z - Target.position.z;
    double dP = sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
    
    //Orientation Error
    double dR = GetAngleDifference(
      EEPose.orientation.x,
      EEPose.orientation.y,
      EEPose.orientation.z,
      EEPose.orientation.w,
      Target.orientation.x,
      Target.orientation.y,
      Target.orientation.z,
      Target.orientation.w);

    //Multi-Objective Weight Randomization
	double random = ZeroOneUniform(RNG);
	
    diffX = EEPose.position.x - BasePose.position.x;
    diffY = EEPose.position.y - BasePose.position.y;
    diffZ = EEPose.position.z - BasePose.position.z;
    double angularScale = sqrt(ChainLength*sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ)) / M_PI;

    return random*dP/angularScale + (1.0-random)*dR;
}

double Evolution::ComputeBalancedFitness(double* &genes) {
	ComputeFK(genes);
	
    //Position Error
    double diffX = EEPose.position.x - Target.position.x;
    double diffY = EEPose.position.y - Target.position.y;
    double diffZ = EEPose.position.z - Target.position.z;
    double dP = sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
    
    //Orientation Error
    double dR = GetAngleDifference(
      EEPose.orientation.x,
      EEPose.orientation.y,
      EEPose.orientation.z,
      EEPose.orientation.w,
      Target.orientation.x,
      Target.orientation.y,
      Target.orientation.z,
      Target.orientation.w);
	
    diffX = EEPose.position.x - BasePose.position.x;
    diffY = EEPose.position.y - BasePose.position.y;
    diffZ = EEPose.position.z - BasePose.position.z;
    double angularScale = sqrt(ChainLength*sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ)) / M_PI;

    return 0.5*(dP/angularScale + dR);
}

double Evolution::ComputeFitness(KDL::Frame &frame) {
	tf::poseKDLToMsg(frame, EEPose);
	
    //Position Error
    double diffX = EEPose.position.x - Target.position.x;
    double diffY = EEPose.position.y - Target.position.y;
    double diffZ = EEPose.position.z - Target.position.z;
    double dP = sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
    
    //Orientation Error
    double dR = GetAngleDifference(
      EEPose.orientation.x,
      EEPose.orientation.y,
      EEPose.orientation.z,
      EEPose.orientation.w,
      Target.orientation.x,
      Target.orientation.y,
      Target.orientation.z,
      Target.orientation.w);

    //Multi-Objective Weight Randomization
	double random = ZeroOneUniform(RNG);
	
    diffX = EEPose.position.x - BasePose.position.x;
    diffY = EEPose.position.y - BasePose.position.y;
    diffZ = EEPose.position.z - BasePose.position.z;
    double angularScale = sqrt(ChainLength*sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ)) / M_PI;

    return random*dP/angularScale + (1.0-random)*dR;
}

void Evolution::ComputeFK(double* &values) {
	KDL::Frame frame;
	int j=0;
    for(unsigned int i=0; i<SegmentCount; i++) {
		if(Chain.getSegment(i).getJoint().getType() != KDL::Joint::None) {
			frame = frame*Chain.getSegment(i).pose(values[j]);
			j++;
		} else {
          frame = frame*Chain.getSegment(i).pose(0.0);
        }
	}
	tf::poseKDLToMsg(frame, EEPose);
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
	double max = Population[PopulationSize-1].Fitness;
	for(int i=0; i<PopulationSize; i++) {
		double grading = (double)i/((double)PopulationSize-1);
		Population[i].Extinction = (Population[i].Fitness + min*(grading-1.0)) / max;
	}
}

bool Evolution::CheckWipeout() {
	//In this, the balanced fitness function with equal weights must be used!
	if(ComputeFitness(Prototype.Genes) < ComputeFitness(Population[0].Genes)) {
		double* values = new double[Dimensionality];
		for(int i=0; i<Dimensionality; i++) {
			values[i] = Population[0].Genes[i];
		}
		for(int i=0; i<Dimensionality; i++) {
			double fitness = ComputeFitness(values);

			double inc = Clip(values[i] + GetRandomValue(0.0, fitness), Dimensions[i]) - values[i];
			values[i] += inc;
			double incFitness = ComputeFitness(values);
			values[i] -= inc;

			double dec = Clip(values[i] - GetRandomValue(0.0, fitness), Dimensions[i]) - values[i];
			values[i] += dec;
			double decFitness = ComputeFitness(values);
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
	double candidateFitness = ComputeBalancedFitness(candidate.Genes);
	if(candidateFitness < EvolutionFitness) {
		EvolutionFitness = candidateFitness;
		for(int i=0; i<Dimensionality; i++) {
			Prototype.Genes[i] = candidate.Genes[i];
			Prototype.Gradient[i] = candidate.Gradient[i];
			Prototype.Extinction = candidate.Extinction;
			Prototype.Fitness = candidate.Fitness;
		}
	}
}

void Evolution::SortByFitness() {
	sort(&Population[0], &Population[PopulationSize], SurvivalOfTheFittest());
}

double Evolution::GetMutationProbability(Individual &parentA, Individual &parentB) {
	double extinction = 0.5 * (parentA.Extinction + parentB.Extinction);
	double inverse = 1.0/Dimensionality;
	return extinction * (1.0-inverse) + inverse;
}

double Evolution::GetMutationStrength(Individual &parentA, Individual &parentB, const Dimension &dimension) {
	double extinction = 0.5 * (parentA.Extinction + parentB.Extinction);
	double span = (dimension.Max - dimension.Min);
	return span * extinction;
}

double Evolution::GetRandomValue(double min, double max) {
	if(min == max) {
		return min;
	}
	boost::random::uniform_real_distribution<> gen(min, max);
	return gen(RNG);
}

int Evolution::GetRandomWeightedIndex(double* &probabilities, int size) {
	double weightSum = 0.0;
	for(int i=0; i<size; i++) {
		weightSum += probabilities[i];
	}
	double rVal = ZeroOneUniform(RNG)*weightSum;
	for(int i=0; i<size; i++) {
		rVal -= probabilities[i];
		if(rVal <= 0.0) {
			return i;
		}
	}
	return size-1;
}

double Evolution::GetAngleDifference(double& q1x, double& q1y, double& q1z, double& q1w, double q2x, double q2y, double q2z, double q2w) {
    if(q1x == q2x && q1y == q2y && q1z == q2z && q1w == q2w) {
      return 0.0;
    }
    double dot = abs(q1x*q2x + q1y*q2y + q1z*q2z + q1w*q2w);
    if(dot > 1.0) {
      dot = 1.0;
    }
    return 2.0 * acos(dot);
}
