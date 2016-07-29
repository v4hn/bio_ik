#include <iostream>
#include <time.h>
#include "ros/ros.h"
#include <tf2/LinearMath/Quaternion.h>

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

  int seed = time(NULL);
  boost::random::mt19937 rng(seed);

  class BIO_IK : public kinematics::KinematicsBase {
    public:
      BIO_IK() {}
      ~BIO_IK() {}

      vector<string> JointNames;
      vector<string> LinkNames;
      vector<string> EndEffector;

      KDL::Chain Chain;
      KDL::JntArray limitsMin, limitsMax;
      int JointCount;

      boost::shared_ptr<KDL::ChainFkSolverPos_recursive> FK_Solver;

      const vector<string>& getJointNames() const { return JointNames; }
      const vector<string>& getLinkNames() const { return LinkNames; }

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

      double IK_FitnessFunction(double* &input, int dimensionality, const geometry_msgs::Pose&) const;
      double IK_BalancedFitnessFunction(const vector<double>& input, int dimensionality, const geometry_msgs::Pose&) const;
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

  double BIO_IK::IK_FitnessFunction(double* &input, int dimensionality, const geometry_msgs::Pose& target) const {
    vector<double> Values(input, input + dimensionality);
    vector<geometry_msgs::Pose> poses;
    getPositionFK(EndEffector, Values, poses);
    
    //Position Error
    double diffX = poses[0].position.x - target.position.x;
    double diffY = poses[0].position.y - target.position.y;
    double diffZ = poses[0].position.z - target.position.z;
    double dP = sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
    
    //Orientation Error
    double dR = GetAngleDifference(
      poses[0].orientation.x,
      poses[0].orientation.y,
      poses[0].orientation.z,
      poses[0].orientation.w,
      target.orientation.x,
      target.orientation.y,
      target.orientation.z,
      target.orientation.w);

    //Multi-Objective Weight Randomization
    boost::random::uniform_real_distribution<> gen(0.0, 1.0);
    double random = gen(rng);

    return random*dP + (1.0-random)*dR;
  }

  double BIO_IK::IK_BalancedFitnessFunction(const vector<double>& input, int dimensionality, const geometry_msgs::Pose& target) const {
    vector<geometry_msgs::Pose> poses;
    getPositionFK(EndEffector, input, poses);
    
    //Position Error
    double diffX = poses[0].position.x - target.position.x;
    double diffY = poses[0].position.y - target.position.y;
    double diffZ = poses[0].position.z - target.position.z;
    double dP = sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
    
    //Orientation Error
    double dR = GetAngleDifference(
      poses[0].orientation.x,
      poses[0].orientation.y,
      poses[0].orientation.z,
      poses[0].orientation.w,
      target.orientation.x,
      target.orientation.y,
      target.orientation.z,
      target.orientation.w);

    //Multi-Objective Weight Randomization
    return 0.5*(dP+dR);
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
    if(!tree.getChain(base_name, tip_name, Chain)) {
      ROS_FATAL("Couldn't find chain %s to %s",base_name.c_str(),tip_name.c_str());
      return false;
    }

    JointCount = Chain.getNrOfJoints();

    limitsMin.resize(JointCount);
    limitsMax.resize(JointCount);

    boost::shared_ptr<const urdf::Joint> joint;
    uint jointNum=0;
    for(unsigned int i=0; i<Chain.segments.size(); ++i) {
      joint = model.getJoint(Chain.segments[i].getJoint().getName());
      if(joint->type != urdf::Joint::UNKNOWN && joint->type != urdf::Joint::FIXED) {
          jointNum++;
          float lower, upper;
          int hasLimits;
          LinkNames.push_back(Chain.segments[i].getName());
          JointNames.push_back(joint->name);

          if( joint->type != urdf::Joint::CONTINUOUS ) {
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
            limitsMin(jointNum-1)=std::numeric_limits<float>::min(); //lowest
            limitsMax(jointNum-1)=std::numeric_limits<float>::max();
          }

          ROS_INFO_STREAM("IK Using joint "<<Chain.segments[i].getName()<<" "<<limitsMin(jointNum-1)<<" "<<limitsMax(jointNum-1));
        }
      }

      EndEffector.push_back(LinkNames[LinkNames.size()-1]);
      FK_Solver.reset(new KDL::ChainFkSolverPos_recursive(Chain));

      return true;
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

    for(unsigned int i=0; i<poses.size(); i++) {
      if(FK_Solver->JntToCart(joints, frame) >= 0) {
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

    solution.resize(JointCount);

    int size = 12;
    int elites = 4;
    int dimensionality = JointCount;
    Dimension* dimensions = new Dimension[dimensionality];
    for(int i=0; i<dimensionality; i++) {
      dimensions[i].Min = limitsMin(i);
      dimensions[i].Max = limitsMax(i);
    }
    Evolution evolution(size, elites, dimensionality, dimensions, boost::bind(&BIO_IK::IK_FitnessFunction, this, _1, _2, ik_pose), ik_seed_state);

    double seedFitness = IK_BalancedFitnessFunction(ik_seed_state, dimensionality, ik_pose);

    double accuracy = 0.001;
    int generations = 0;
    while((double)(clock() - begin_time) / CLOCKS_PER_SEC < 0.02) {
      evolution.Evolve();
      generations += 1;
    }
    for(int i=0; i<JointCount; i++) {
      solution[i] = evolution.GetPrototype().Genes[i];
    }

    double solutionFitness = IK_BalancedFitnessFunction(solution, dimensionality, ik_pose);

    if(seedFitness <= solutionFitness) {
      solution = ik_seed_state;
      cout << "Accuracy: " << seedFitness << endl;
    } else {
      cout << "Accuracy: " << solutionFitness << endl;
    }
    //cout << "Generations: " << generations << " (" << ((double)(clock() - begin_time) / CLOCKS_PER_SEC) << "s)" << endl;
    return true;
  }

}

//register BIO_IK as a KinematicsBase implementation
#include <pluginlib/class_list_macros.h>
PLUGINLIB_EXPORT_CLASS(bio_ik::BIO_IK, kinematics::KinematicsBase);
