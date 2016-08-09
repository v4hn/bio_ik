#include <iostream>
#include "bio_ik.cpp"

double GetAverage(vector<double> &values) {
  double sum = 0.0;
  for(int i=0; i<values.size(); i++) {
    sum += values[i];
  }
  return sum / values.size();
}

double GetMaximum(vector<double> &values) {
  double max = 0.0;
  for(int i=0; i<values.size(); i++) {
    if(values[i] > max) {
      max = values[i];
    }
  }
  return max;
}

double GetMinimum(vector<double> &values) {
  double min = 1000.0;
  for(int i=0; i<values.size(); i++) {
    if(values[i] < min) {
      min = values[i];
    }
  }
  return min;
}

void FilterSolutions(vector< vector<double> > solutions, vector<bool> success, vector< vector<double> > &filtered, double threshold) {
  for(int i=0; i<solutions.size(); i++) {
    if(success[i]) {
      if(filtered.size() > 0) {
        bool add = true;
        for(int j=0; j<filtered.size(); j++) {
          double dist = 0.0;
          for(int k=0; k<filtered[j].size(); k++) {
            double diff = filtered[j][k] - solutions[i][k];
            dist += diff*diff;
          }
          dist = sqrt(dist);
          if(dist < threshold) {
            add = false;
          }
        }
        if(add) {
          filtered.push_back(solutions[i]);
        }
      } else {
        filtered.push_back(solutions[i]);
      }
    }
  }
}

void FilterOperation(bio_ik::BIO_IK &solver) {
  vector< vector<double> > Solutions;
  vector<double> FitnessResults;
  vector<bool> Success;
  vector<double> ComputationTimes;
  vector<double> GenerationCounts;
  
  /*
  vector<string> links;
  //links.push_back("r_wrist_roll_link");
  links.push_back("drill_chunk_drill_mount");

  vector<double> values;
  solver.GetRandomConfiguration(values);

  //vector<geometry_msgs::Pose> poses;
  //solver.getPositionFK(links, values, poses);
  geometry_msgs::Pose pose = poses[0];
  */
  /*
  geometry_msgs::Pose pose;

  pose.position.x = 0.0689272;
  pose.position.y = -0.045;
  pose.position.z = 1.422;
  pose.orientation.x = -0.000;
  pose.orientation.y = 0.980;
  pose.orientation.z = 0.003;
  pose.orientation.w = 0.198;
  */
  vector<string> links;
  links.push_back("pa10/pa10_T6_link");
  vector<double> values;
  values.push_back(2.6517);
  values.push_back(-0.8852);
  values.push_back(1.3951);
  values.push_back(-2.7008);
  values.push_back(1.6079);
  values.push_back(0.2346);

  //vector<geometry_msgs::Pose> poses;
  vector<geometry_msgs::Pose> poses;
  solver.getPositionFK(links, values, poses);
  geometry_msgs::Pose pose = poses[0];

  for(int i=0; i<100; i++) {
    vector<double> seed;
    //seed.resize(solver.JointCount);
    solver.GetRandomConfiguration(seed);
    vector<double> solution;
    double solutionFitness;
    double computationTime;
    double generationCount;

    Success.push_back(solver.myPositionIK(pose, seed, solution, solutionFitness, computationTime, generationCount, 0.001, 0.1));
    FitnessResults.push_back(solutionFitness);
    ComputationTimes.push_back(computationTime);
    GenerationCounts.push_back(generationCount);
    Solutions.push_back(solution);
  }

  double avgFitness = GetAverage(FitnessResults);
  double successRate = 0.0;
  for(int i=0; i<Success.size(); i++) {
    if(Success[i]) {
      successRate += 1.0;
    }
  }
  successRate /= Success.size();
  double avgTime = GetAverage(ComputationTimes);
  double minTime = GetMinimum(ComputationTimes);
  double maxTime = GetMaximum(ComputationTimes);

  cout << "Average Fitness: " << avgFitness << endl;
  cout << "Success Rate: " << successRate << endl;
  cout << "Average Time: " << avgTime << endl;
  cout << "Min Time: " << minTime << endl;
  cout << "Max Time: " << maxTime << endl;

  cout << "ALL: " << Solutions.size() << endl;
  for(int i=0; i<Solutions.size() ;i++) {
    cout << "Solution " << i << " : ";
    for(int j=0; j<Solutions[i].size(); j++) {
      cout << Solutions[i][j] << " ";
    }
    cout << endl;
  }

  vector< vector<double> > filtered;
  FilterSolutions(Solutions, Success, filtered, 2.5);
  cout << "FILTERED: " << filtered.size() << endl;

  for(int i=0; i<filtered.size() ;i++) {
    cout << "Filtered " << i << " : ";
    for(int j=0; j<filtered[i].size(); j++) {
      cout << filtered[i][j] << " ";
    }
    cout << endl;
  }
}

void SearchOperation(bio_ik::BIO_IK solver, int sampleCount) {
  vector< vector<double> > Solutions;
  vector<double> FitnessResults;
  vector<bool> Success;
  vector<double> ComputationTimes;
  vector<double> GenerationCounts;
  
  vector<string> links;
  links.push_back("r_wrist_roll_link");
  //links.push_back("drill_chunk_drill_mount");

  for(int k=0; k<sampleCount; k++) {
    vector<double> values;
    solver.GetRandomConfiguration(values);

    vector<geometry_msgs::Pose> poses;
    solver.getPositionFK(links, values, poses);
    geometry_msgs::Pose pose = poses[0];

    vector<double> seed;
    seed.resize(solver.JointCount);
    for(int i=0; i<seed.size(); i++) {
      seed[i] = 0.0;
    }
    //solver.GetRandomConfiguration(seed);
    vector<double> solution;
    double solutionFitness;
    double computationTime;
    double generationCount;

    Success.push_back(solver.myPositionIK(pose, seed, solution, solutionFitness, computationTime, generationCount, 0.001, 10.0));
    FitnessResults.push_back(solutionFitness);
    ComputationTimes.push_back(computationTime);
    GenerationCounts.push_back(generationCount);
    Solutions.push_back(solution);
    
    cout << "Optimized sample " << k << "!" << endl;
  }

  double avgFitness = GetAverage(FitnessResults);
  double successRate = 0.0;
  for(int i=0; i<Success.size(); i++) {
    if(Success[i]) {
      successRate += 1.0;
    }
  }
  successRate /= Success.size();
  double avgTime = GetAverage(ComputationTimes);
  double minTime = GetMinimum(ComputationTimes);
  double maxTime = GetMaximum(ComputationTimes);
  double avgGenerations = GetAverage(GenerationCounts);

  cout << "=========================" << endl,
  cout << "Samples: " << sampleCount << endl;
  cout << "Average Fitness: " << avgFitness << endl;
  cout << "Success Rate: " << successRate << endl;
  cout << "Average Time: " << avgTime << endl;
  cout << "Min Time: " << minTime << endl;
  cout << "Max Time: " << maxTime << endl;
  cout << "Average Generations: " << avgGenerations << endl;
  cout << "=========================" << endl;
}

int main(int argc, char *argv[]) {
  ros::init(argc, argv, "bio_ik", 1);
  
  bio_ik::BIO_IK solver;
  solver.initialize("robot_description", "right_arm", "torso_lift_link", "r_wrist_roll_link", 0.0);
  //solver.initialize("robot_description", "pa10_planning_group", "world", "pa10/pa10_T6_link", 0.0);

  //FilterOperation(solver);
  SearchOperation(solver, 1000);

  return 0;
}
