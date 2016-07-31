#include <iostream>
#include "bio_ik.cpp"

int main(int argc, char *argv[]) {
  ros::init(argc, argv, "bio_ik");
  
  bio_ik::BIO_IK solver;
  solver.initialize("robot_description", "right_arm", "base_link", "r_wrist_roll_link", 0.0);

  vector<string> links;
  links.push_back("r_wrist_roll_link");
  vector<double> values;
  solver.GetRandomConfiguration(values);
  
  vector<geometry_msgs::Pose> poses;
  solver.getPositionFK(links, values, poses);
  geometry_msgs::Pose pose = poses[0];

  vector<double> seed;
  seed.resize(solver.JointCount);
  vector<double> solution;
  moveit_msgs::MoveItErrorCodes err;

  solver.getPositionIK(pose, seed, solution, err);

  //122-148 Generations at 0.01s

  return 0;
}
