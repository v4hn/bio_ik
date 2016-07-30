#include <iostream>
#include <time.h>
#include "ros/ros.h"

#include "evolution.h"

#include <moveit/kinematics_base/kinematics_base.h>

#include <kdl_parser/kdl_parser.hpp>
#include <kdl/chain.hpp>
#include <kdl/tree.hpp>
#include "kdl/frames.hpp"
#include "kdl/jntarray.hpp"

#include "kdl/chainfksolverpos_recursive.hpp"

#include <tf_conversions/tf_kdl.h>

#include <urdf/model.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

using namespace std;

namespace bio_ik {

  boost::random::mt19937 rng(time(NULL));

  class BIO_IK : public kinematics::KinematicsBase {
    public:
      BIO_IK() {}
      ~BIO_IK() {}

      vector<string> JointNames, LinkNames;

      int IndexBase, IndexEE;

      KDL::Chain Chain;
      double ChainLength;
      KDL::JntArray limitsMin, limitsMax;
      int JointCount;
      int SegmentCount;

      const vector<string>& getJointNames() const { return JointNames; }
      const vector<string>& getLinkNames() const { return LinkNames; }

      void getPositionFK_BioIK(const vector<double> &joint_angles,
                              geometry_msgs::Pose &eePose,
                              geometry_msgs::Pose &basePose) const;

      void getPositionFK_BioIK(double* &joint_angles,
                              geometry_msgs::Pose &eePose,
                              geometry_msgs::Pose &basePose) const;

      bool getPositionFK(const vector<string> &link_names,
                        const vector<double> &joint_angles,
                        vector<geometry_msgs::Pose> &poses) const;

      bool getPositionIK(const geometry_msgs::Pose &ik_pose,
                        const vector<double> &ik_seed_state,
                        vector<double> &solution,
                        moveit_msgs::MoveItErrorCodes &error_code,
                        const kinematics::KinematicsQueryOptions &options = kinematics::KinematicsQueryOptions()) const;

      bool searchPositionIK(const geometry_msgs::Pose &ik_pose,
                            const vector<double> &ik_seed_state,
                            double timeout,
                            vector<double> &solution,
                            moveit_msgs::MoveItErrorCodes &error_code,
                            const kinematics::KinematicsQueryOptions &options = kinematics::KinematicsQueryOptions()) const;

      bool searchPositionIK(const geometry_msgs::Pose &ik_pose,
                            const vector<double> &ik_seed_state,
                            double timeout,
                            const vector<double> &consistency_limits,
                            vector<double> &solution,
                            moveit_msgs::MoveItErrorCodes &error_code,
                            const kinematics::KinematicsQueryOptions &options = kinematics::KinematicsQueryOptions()) const;

      bool searchPositionIK(const geometry_msgs::Pose &ik_pose,
                            const vector<double> &ik_seed_state,
                            double timeout,
                            vector<double> &solution,
                            const IKCallbackFn &solution_callback,
                            moveit_msgs::MoveItErrorCodes &error_code,
                            const kinematics::KinematicsQueryOptions &options = kinematics::KinematicsQueryOptions()) const;

      bool searchPositionIK(const geometry_msgs::Pose &ik_pose,
                            const vector<double> &ik_seed_state,
                            double timeout,
                            const vector<double> &consistency_limits,
                            vector<double> &solution,
                            const IKCallbackFn &solution_callback,
                            moveit_msgs::MoveItErrorCodes &error_code,
                            const kinematics::KinematicsQueryOptions &options = kinematics::KinematicsQueryOptions()) const;

      bool searchPositionIK(const geometry_msgs::Pose &ik_pose,
                            const vector<double> &ik_seed_state,
                            double timeout,
                            vector<double> &solution,
                            const IKCallbackFn &solution_callback,
                            moveit_msgs::MoveItErrorCodes &error_code,
                            const vector<double> &consistency_limits,
                            const kinematics::KinematicsQueryOptions &options) const;

      bool initialize(const string &robot_description,
                      const string& group_name,
                      const string& base_name,
                      const string& tip_name,
                      double search_discretization);

      int getKDLSegmentIndex(const string &name) const;

      void GetRandomConfiguration(vector<double>& configuration);

      double IK_FitnessFunction(double* &input, const geometry_msgs::Pose&) const;
      double IK_BalancedFitnessFunction(const vector<double>& input, const geometry_msgs::Pose&) const;
    private:

  };

  double GetAngleDifference(double& q1x, double& q1y, double& q1z, double& q1w, double q2x, double q2y, double q2z, double q2w) {
    if(q1x == q2x && q1y == q2y && q1z == q2z && q1w == q2w) {
      return 0.0;
    }
    double dot = abs(q1x*q2x + q1y*q2y + q1z*q2z + q1w*q2w);
    if(dot > 1.0) {
      dot = 1.0;
    }
    return 2.0 * acos(dot);
  }

  double BIO_IK::IK_FitnessFunction(double* &input, const geometry_msgs::Pose& target) const {
    geometry_msgs::Pose eePose, basePose;
    getPositionFK_BioIK(input, eePose, basePose);
    
    //Position Error
    double diffX = eePose.position.x - target.position.x;
    double diffY = eePose.position.y - target.position.y;
    double diffZ = eePose.position.z - target.position.z;
    double dP = sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
    
    //Orientation Error
    double dR = GetAngleDifference(
      eePose.orientation.x,
      eePose.orientation.y,
      eePose.orientation.z,
      eePose.orientation.w,
      target.orientation.x,
      target.orientation.y,
      target.orientation.z,
      target.orientation.w);

    //Multi-Objective Weight Randomization
    boost::random::uniform_real_distribution<> gen(0.0, 1.0);
    double random = gen(rng);

    diffX = eePose.position.x - basePose.position.x;
    diffY = eePose.position.y - basePose.position.y;
    diffZ = eePose.position.z - basePose.position.z;
    double angularScale = sqrt(ChainLength*sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ)) / M_PI;

    return random*dP/angularScale + (1.0-random)*dR;
  }

  double BIO_IK::IK_BalancedFitnessFunction(const vector<double>& input, const geometry_msgs::Pose& target) const {
    geometry_msgs::Pose eePose, basePose;
    getPositionFK_BioIK(input, eePose, basePose);
    
    //Position Error
    double diffX = eePose.position.x - target.position.x;
    double diffY = eePose.position.y - target.position.y;
    double diffZ = eePose.position.z - target.position.z;
    double dP = sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
    
    //Orientation Error
    double dR = GetAngleDifference(
      eePose.orientation.x,
      eePose.orientation.y,
      eePose.orientation.z,
      eePose.orientation.w,
      target.orientation.x,
      target.orientation.y,
      target.orientation.z,
      target.orientation.w);

    //Multi-Objective Weight Randomization
    diffX = eePose.position.x - basePose.position.x;
    diffY = eePose.position.y - basePose.position.y;
    diffZ = eePose.position.z - basePose.position.z;
    double angularScale = sqrt(ChainLength*sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ)) / M_PI;
    
    return 0.5*(dP/angularScale + dR);
  }

  bool BIO_IK::initialize(const string &robot_description,
                          const string& group_name, 
                          const string& base_name, 
                          const string& tip_name, 
                          double search_discretization) {
    setValues(robot_description, group_name, base_name, tip_name, search_discretization);
    ros::NodeHandle nh("~");

    urdf::Model model;
    string xml_string, urdf_xml, full_urdf_xml;

    nh.param("urdf_xml",urdf_xml,robot_description);
    nh.searchParam(urdf_xml,full_urdf_xml);
    
    ROS_DEBUG_NAMED("bio_ik","Reading xml file from parameter server");
    if (!nh.getParam(full_urdf_xml, xml_string)) {
      ROS_FATAL_NAMED("bio_ik","Could not load the xml from parameter server: %s", urdf_xml.c_str());
      return false;
    }
    nh.param(full_urdf_xml,xml_string,string());
    model.initString(xml_string);

    ROS_DEBUG_STREAM_NAMED("bio_ik","Reading joints and links from URDF");
    KDL::Tree tree;
    if(!kdl_parser::treeFromUrdfModel(model, tree)) {
      ROS_FATAL("Failed to extract kdl tree from xml robot description");
      return false;
    }
    cout << "Using chain from " << base_name << " to " << tip_name << "!" << endl;
    if(!tree.getChain(base_name, tip_name, Chain)) {
      ROS_FATAL("Couldn't find chain %s to %s",base_name.c_str(),tip_name.c_str());
      return false;
    }

    SegmentCount = Chain.getNrOfSegments();
    JointCount = Chain.getNrOfJoints();

    limitsMin.resize(JointCount);
    limitsMax.resize(JointCount);

    boost::shared_ptr<const urdf::Joint> joint;
    uint jointNum=0;
    for(unsigned int i=0; i<SegmentCount; ++i) {
      joint = model.getJoint(Chain.segments[i].getJoint().getName());
      if(joint->type != urdf::Joint::UNKNOWN && joint->type != urdf::Joint::FIXED) {
          jointNum++;
          float lower, upper;
          int hasLimits;
          LinkNames.push_back(Chain.segments[i].getName());
          JointNames.push_back(joint->name);
          ChainLength += Chain.segments[i].getJoint().JointOrigin().Norm();

          if(joint->type != urdf::Joint::CONTINUOUS) {
            if(joint->safety) {
              lower = std::max(joint->limits->lower, joint->safety->soft_lower_limit);
              upper = std::min(joint->limits->upper, joint->safety->soft_upper_limit);
            } else {
              lower = joint->limits->lower;
              upper = joint->limits->upper;
            }
            hasLimits = 1;
          }
          else {
            hasLimits = 0;
          }

          if(hasLimits) {
            limitsMin(jointNum-1)=lower;
            limitsMax(jointNum-1)=upper;
          }
          else {
            limitsMin(jointNum-1)=std::numeric_limits<float>::min();
            limitsMax(jointNum-1)=std::numeric_limits<float>::max();
          }

          ROS_INFO_STREAM("IK Using joint "<<Chain.segments[i].getName()<<" "<<limitsMin(jointNum-1)<<" "<<limitsMax(jointNum-1));
        }
      }

      IndexBase = getKDLSegmentIndex(LinkNames[0]);
      IndexEE = getKDLSegmentIndex(LinkNames[LinkNames.size()-1]);
      
      return true;
  }

  int BIO_IK::getKDLSegmentIndex(const string &name) const {
    int i=0;
    while(i < SegmentCount) {
      if(Chain.getSegment(i).getName() == name) {
        return i+1;
      }
      i++;
    }
    return -1;
  }

  void BIO_IK::GetRandomConfiguration(vector<double>& configuration) {
    configuration.resize(JointCount);
    for(int i=0; i<JointCount; i++) {
      boost::random::uniform_real_distribution<> gen(limitsMin(i), limitsMax(i));
      configuration[i] = gen(rng);
    }
  }

  void BIO_IK::getPositionFK_BioIK(const vector<double> &joint_angles,
                                  geometry_msgs::Pose &eePose,
                                  geometry_msgs::Pose &basePose) const {
      KDL::Frame frame;
      int j=0;
      for(unsigned int i=0; i<SegmentCount; i++) {
        if(Chain.getSegment(i).getJoint().getType() != KDL::Joint::None) {
          frame = frame*Chain.getSegment(i).pose(joint_angles[j]);
          j++;
          if(j == IndexBase) {
            tf::poseKDLToMsg(frame, basePose);
          }
        } else {
          frame = frame*Chain.getSegment(i).pose(0.0);
        }
      }
      tf::poseKDLToMsg(frame, eePose);
  }

  void BIO_IK::getPositionFK_BioIK(double* &joint_angles,
                                  geometry_msgs::Pose &eePose,
                                  geometry_msgs::Pose &basePose) const {
      KDL::Frame frame;
      int j=0;
      for(unsigned int i=0; i<SegmentCount; i++) {
        if(Chain.getSegment(i).getJoint().getType() != KDL::Joint::None) {
          frame = frame*Chain.getSegment(i).pose(joint_angles[j]);
          j++;
          if(j == IndexBase) {
            tf::poseKDLToMsg(frame, basePose);
          }
        } else {
          frame = frame*Chain.getSegment(i).pose(0.0);
        }
      }
      tf::poseKDLToMsg(frame, eePose);
  }

  bool BIO_IK::getPositionFK(const vector<string> &link_names, 
                            const vector<double> &joint_angles, 
                            vector<geometry_msgs::Pose> &poses) const {
    poses.resize(link_names.size());
    if(joint_angles.size() != JointCount) {
      cout << "Joint angle size incorrect." << endl;
      return false;
    }

    KDL::JntArray joints(JointCount);
    KDL::Frame frame;
    for(unsigned int i=0; i<JointCount; i++) {
      joints(i) = joint_angles[i];
    }

    KDL::ChainFkSolverPos_recursive fk_solver(Chain);
    for(unsigned int i=0; i<poses.size(); i++) {
      if(fk_solver.JntToCart(joints, frame, getKDLSegmentIndex(link_names[i])) >= 0) {
        tf::poseKDLToMsg(frame, poses[i]);
      } else {
        cout << "Failed to compute FK for joint " << i << endl;
        return false;
      }
    }

    return true;
  }

  bool BIO_IK::getPositionIK(const geometry_msgs::Pose &ik_pose,
                            const vector<double> &ik_seed_state,
                            vector<double> &solution,
                            moveit_msgs::MoveItErrorCodes &error_code,
                            const kinematics::KinematicsQueryOptions &options) const {
    const IKCallbackFn solution_callback = 0;
    std::vector<double> consistency_limits;

    return searchPositionIK(ik_pose,
                            ik_seed_state,
                            default_timeout_,
                            solution,
                            solution_callback,
                            error_code,
                            consistency_limits,
                            options);
  }

  bool BIO_IK::searchPositionIK(const geometry_msgs::Pose &ik_pose,
                                const vector<double> &ik_seed_state,
                                double timeout,
                                vector<double> &solution,
                                moveit_msgs::MoveItErrorCodes &error_code,
                                const kinematics::KinematicsQueryOptions &options) const {
    const IKCallbackFn solution_callback = 0;
    std::vector<double> consistency_limits;

    return searchPositionIK(ik_pose,
                            ik_seed_state,
                            timeout,
                            solution,
                            solution_callback,
                            error_code,
                            consistency_limits,
                            options);
  }

  bool BIO_IK::searchPositionIK(const geometry_msgs::Pose &ik_pose,
                                const vector<double> &ik_seed_state,
                                double timeout,
                                const vector<double> &consistency_limits,
                                vector<double> &solution,
                                moveit_msgs::MoveItErrorCodes &error_code,
                                const kinematics::KinematicsQueryOptions &options) const {
    const IKCallbackFn solution_callback = 0;
    return searchPositionIK(ik_pose,
                            ik_seed_state,
                            timeout,
                            solution,
                            solution_callback,
                            error_code,
                            consistency_limits,
                            options);
  }

  bool BIO_IK::searchPositionIK(const geometry_msgs::Pose &ik_pose,
                                const vector<double> &ik_seed_state,
                                double timeout,
                                vector<double> &solution,
                                const IKCallbackFn &solution_callback,
                                moveit_msgs::MoveItErrorCodes &error_code,
                                const kinematics::KinematicsQueryOptions &options) const {
    std::vector<double> consistency_limits;
    return searchPositionIK(ik_pose,
                            ik_seed_state,
                            timeout,
                            solution,
                            solution_callback,
                            error_code,
                            consistency_limits,
                            options);
  }

  bool BIO_IK::searchPositionIK(const geometry_msgs::Pose &ik_pose,
                                const vector<double> &ik_seed_state,
                                double timeout,
                                const vector<double> &consistency_limits,
                                vector<double> &solution,
                                const IKCallbackFn &solution_callback,
                                moveit_msgs::MoveItErrorCodes &error_code,
                                const kinematics::KinematicsQueryOptions &options) const {
    return searchPositionIK(ik_pose,
                            ik_seed_state,
                            timeout,
                            solution,
                            solution_callback,
                            error_code,
                            consistency_limits,
                            options);
  }

  bool BIO_IK::searchPositionIK(const geometry_msgs::Pose &ik_pose,
                                const vector<double> &ik_seed_state,
                                double timeout,
                                vector<double> &solution,
                                const IKCallbackFn &solution_callback,
                                moveit_msgs::MoveItErrorCodes &error_code,
                                const vector<double> &consistency_limits,
                                const kinematics::KinematicsQueryOptions &options) const {
    const clock_t begin_time = clock();

    ROS_DEBUG_STREAM_NAMED("bio_ik","getPositionIK");

    KDL::Frame frame;
    tf::poseMsgToKDL(ik_pose,frame);

    Dimension* dimensions = new Dimension[JointCount];
    for(int i=0; i<JointCount; i++) {
      dimensions[i].Min = limitsMin(i);
      dimensions[i].Max = limitsMax(i);
    }
    Evolution evolution(12, 3, JointCount, dimensions, boost::bind(&BIO_IK::IK_FitnessFunction, this, _1, ik_pose), ik_seed_state);
    double seedFitness = IK_BalancedFitnessFunction(ik_seed_state, ik_pose);

    //double accuracy = 0.001;
    int generations = 0;
    while((double)(clock() - begin_time) / CLOCKS_PER_SEC < 0.01) {
      evolution.Evolve();
      generations += 1;
    }

    solution.resize(JointCount);
    for(int i=0; i<JointCount; i++) {
      solution[i] = evolution.GetPrototype().Genes[i];
    }
    double solutionFitness = IK_BalancedFitnessFunction(solution, ik_pose);

    if(seedFitness <= solutionFitness) {
        solution = ik_seed_state;
        cout << "Generations: " << generations << " Fitness: " <<  seedFitness << " (" << ((double)(clock() - begin_time) / CLOCKS_PER_SEC) << "s)" << endl;
    } else {
        cout << "Generations: " << generations << " Fitness: " <<  solutionFitness << " (" << ((double)(clock() - begin_time) / CLOCKS_PER_SEC) << "s)" << endl;
    }
    
    return true;
  }

}

//register BIO_IK as a KinematicsBase implementation
#include <pluginlib/class_list_macros.h>
PLUGINLIB_EXPORT_CLASS(bio_ik::BIO_IK, kinematics::KinematicsBase);
