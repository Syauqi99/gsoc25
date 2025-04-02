cmake_minimum_required(VERSION 3.10)
project(GeomScaleGSoC25)

set(CMAKE_CXX_STANDARD 14)

# Locate dependencies
find_package(Eigen3 REQUIRED)
find_package(Volesti REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR} ${Volesti_INCLUDE_DIRS})

# Download Volesti package into third_partylibrary
include(ExternalProject)
ExternalProject_Add(volesti
    GIT_REPOSITORY https://github.com/volesti/volesti.git
    GIT_TAG master
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_partylibrary/volesti
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
