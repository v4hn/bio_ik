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

      int JointCount, SegmentCount;
      int IndexBase, IndexEE;
      KDL::Chain Chain;
      double ChainLength;
      Dimension* Dimensions;

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

      int getKDLSegmentIndex(const string &name) const;

      void GetRandomConfiguration(vector<double>& configuration);

    private:

  };

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
    Dimensions = new Dimension[JointCount];

    boost::shared_ptr<const urdf::Joint> joint;
    uint jointNum=0;
    for(unsigned int i=0; i<SegmentCount; ++i) {
      joint = model.getJoint(Chain.segments[i].getJoint().getName());
      
      if(joint->type != urdf::Joint::UNKNOWN && joint->type != urdf::Joint::FIXED) {
          jointNum++;

          LinkNames.push_back(Chain.segments[i].getName());
          JointNames.push_back(joint->name);

          ChainLength += Chain.segments[i].getJoint().JointOrigin().Norm();

          if(joint->type != urdf::Joint::CONTINUOUS) {
            if(joint->safety) {
              Dimensions[jointNum-1].Min = max(joint->limits->lower, joint->safety->soft_lower_limit);
              Dimensions[jointNum-1].Max = min(joint->limits->upper, joint->safety->soft_upper_limit);
            } else {
              Dimensions[jointNum-1].Min = joint->limits->lower;
              Dimensions[jointNum-1].Max = joint->limits->upper;
            }
          } else {
            Dimensions[jointNum-1].Min = numeric_limits<float>::min();
            Dimensions[jointNum-1].Max = numeric_limits<float>::max();
          }

          ROS_INFO_STREAM("IK Using joint "<<Chain.segments[i].getName()<<" "<<Dimensions[jointNum-1].Min<<" "<<Dimensions[jointNum-1].Max);
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
      boost::random::uniform_real_distribution<> gen(Dimensions[i].Min, Dimensions[i].Max);
      configuration[i] = gen(rng);
    }
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

    Evolution evolution(12, 3, JointCount, Dimensions, ik_seed_state, ik_pose, Chain, ChainLength, JointCount, SegmentCount);
    
    //double accuracy = 0.001;
    int generations = 0;
    while((double)(clock() - begin_time) / CLOCKS_PER_SEC < 0.005) {
      evolution.Evolve();
      generations += 1;
    }

    solution.resize(JointCount);
    for(int i=0; i<JointCount; i++) {
      solution[i] = evolution.GetPrototype().Genes[i];
    }

    cout << "Generations: " << generations << " Evolution Fitness: " <<  evolution.GetEvolutionFitness() << " (" << ((double)(clock() - begin_time) / CLOCKS_PER_SEC) << "s)" << endl;
    
    return true;
  }

}

//register BIO_IK as a KinematicsBase implementation
#include <pluginlib/class_list_macros.h>
PLUGINLIB_EXPORT_CLASS(bio_ik::BIO_IK, kinematics::KinematicsBase);
