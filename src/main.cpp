#include <iostream>
#include "bio_ik.cpp"

void FilterSolutions(vector< vector<double> > solutions, vector< vector<double> > &filtered, double threshold) {
  for(int i=0; i<solutions.size(); i++) {
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

int main(int argc, char *argv[]) {
  ros::init(argc, argv, "bio_ik");
  
  bio_ik::BIO_IK solver;
  //solver.initialize("robot_description", "right_arm", "base_link", "r_wrist_roll_link", 0.0);
  solver.initialize("robot_description", "pa10_planning_group", "pa10/pa10_s1_link", "drill_chunk_drill_mount", 0.0);

  vector< vector<double> > Solutions;
  vector<double> FitnessResults;
  vector<bool> Success;
  vector<double> ComputationTimes;
  double MinTime = 1000.0;
  double MaxTime = 0.0;
  
  vector<string> links;
  //links.push_back("r_wrist_roll_link");
  links.push_back("drill_chunk_drill_mount");

  vector<double> values;
  solver.GetRandomConfiguration(values);

  vector<geometry_msgs::Pose> poses;
  solver.getPositionFK(links, values, poses);
  geometry_msgs::Pose pose = poses[0];

  /*
  geometry_msgs::Pose pose;
  pose.position.x = -0.248;
  pose.position.y = -0.045;
  pose.position.z = 1.422;
  pose.orientation.x = -0.000;
  pose.orientation.y = 0.980;
  pose.orientation.z = 0.003;
  pose.orientation.w = 0.198;
  */

  for(int i=0; i<1000; i++) {
    vector<double> seed;
    //seed.resize(solver.JointCount);
    solver.GetRandomConfiguration(seed);
    vector<double> solution;
    double solutionFitness;
    double computationTime;

    Success.push_back(solver.myPositionIK(pose, seed, solution, solutionFitness, computationTime, 0.001, 0.1));
    FitnessResults.push_back(solutionFitness);
    ComputationTimes.push_back(computationTime);
    Solutions.push_back(solution);
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

  cout << "ALL: " << Solutions.size() << endl;
  for(int i=0; i<Solutions.size() ;i++) {
    cout << "Solution " << i << " : ";
    for(int j=0; j<Solutions[i].size(); j++) {
      cout << Solutions[i][j] << " ";
    }
    cout << endl;
  }

  vector< vector<double> > filtered;
  FilterSolutions(Solutions, filtered, 2.5);
  cout << "FILTERED: " << filtered.size() << endl;

  for(int i=0; i<filtered.size() ;i++) {
    cout << "Filtered " << i << " : ";
    for(int j=0; j<filtered[i].size(); j++) {
      cout << filtered[i][j] << " ";
    }
    cout << endl;
  }

  return 0;
}
