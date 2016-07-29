#include <iostream>
#include "bio_ik.cpp"

int main(int argc, char *argv[]) {
  ros::init(argc, argv, "bio_ik");
  
  bio_ik::BIO_IK solver;
  solver.initialize("robot_description", "right_arm", "base_link", "r_wrist_roll_link", 0.0);

  geometry_msgs::Pose target;
  target.orientation.w = 1.0;
  vector<double> seed;
  for(int i=0; i<solver.Chain.getNrOfJoints(); i++) {
    seed.push_back(0.0);
  }
  vector<double> solution;
  moveit_msgs::MoveItErrorCodes err;
  solver.getPositionIK(target, seed, solution, err);

  return 0;
}
