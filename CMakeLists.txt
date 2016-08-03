cmake_minimum_required(VERSION 2.8.3)
project(bio_ik)

find_package(catkin REQUIRED COMPONENTS
    roscpp
    moveit_core
    pluginlib
    tf_conversions
    kdl_parser
)

catkin_package(
  CATKIN_DEPENDS
    roscpp
    moveit_core
    pluginlib
    tf_conversions
)

include_directories(${catkin_INCLUDE_DIRS})

add_library(bio_ik
    src/bio_ik.cpp
    src/evolution.cpp
)

target_link_libraries(bio_ik
    ${catkin_LIBRARIES}
)

add_executable(bio_ik_test
    src/main.cpp
)

target_link_libraries(bio_ik_test
    ${catkin_LIBRARIES}
    bio_ik
)
