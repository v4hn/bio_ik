#include <iostream>
#include "bio_ik.cpp"

int main(int argc, char *argv[]) {
  ros::init(argc, argv, "bio_ik");
  
  bio_ik::BIO_IK solver;
  solver.initialize("robot_description", "right_arm", "base_link", "r_wrist_roll_link", 0.0);

  vector<double> FitnessResults;
  vector<bool> Success;
  vector<double> ComputationTimes;
  double MinTime = 1000.0;
  double MaxTime = 0.0;

  vector<double> seed;
  seed.resize(solver.JointCount);
  
  vector<string> links;
  links.push_back("r_wrist_roll_link");
  
  for(int i=0; i<1; i++) {
    vector<double> values;
    solver.GetRandomConfiguration(values);
    
    vector<geometry_msgs::Pose> poses;
    solver.getPositionFK(links, values, poses);
    geometry_msgs::Pose pose = poses[0];

    vector<double> solution;

    double solutionFitness;

    double computationTime;

    Success.push_back(solver.myPositionIK(pose, seed, solution, solutionFitness, computationTime));
    FitnessResults.push_back(solutionFitness);
    ComputationTimes.push_back(computationTime);
  }

  double avgFitness = 0.0;
  for(int i=0; i<FitnessResults.size(); i++) {
    avgFitness += FitnessResults[i];
  }
  avgFitness /= FitnessResults.size();

  double successRate = 0.0;
  for(int i=0; i<Success.size(); i++) {
    if(Success[i]) {
      successRate += 1.0;
    }
  }
  successRate /= Success.size();

  double avgTime = 0.0;
  for(int i=0; i<ComputationTimes.size(); i++) {
    avgTime += ComputationTimes[i];
    if(ComputationTimes[i] < MinTime) {
      MinTime = ComputationTimes[i];
    }
    if(ComputationTimes[i] > MaxTime) {
      MaxTime = ComputationTimes[i];
    }
  }
  avgTime /= ComputationTimes.size();

  cout << "Average Fitness: " << avgFitness << endl;
  cout << "Success Rate: " << successRate << endl;
  cout << "Average Time: " << avgTime << endl;
  cout << "Min Time: " << MinTime << endl;
  cout << "Max Time: " << MaxTime << endl;

  return 0;
}
