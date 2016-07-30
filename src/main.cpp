#include <iostream>
#include "bio_ik.cpp"

int main(int argc, char *argv[]) {
  ros::init(argc, argv, "bio_ik");
  
  bio_ik::BIO_IK solver;
  solver.initialize("robot_description", "right_arm", "base_link", "r_wrist_roll_link", 0.0);

  geometry_msgs::Pose pose;
  pose.position.x = 0.638568;
  pose.position.y = -0.487843;
  pose.position.z = 0.867719;
  pose.orientation.x = 0.629086;
  pose.orientation.y = 0.469335;
  pose.orientation.z = -0.390434;
  pose.orientation.w = 0.481183;

  vector<double> seed;
  seed.resize(solver.JointCount);

  vector<double> solution;

  moveit_msgs::MoveItErrorCodes err;

  solver.getPositionIK(pose, seed, solution, err);

  return 0;
}
