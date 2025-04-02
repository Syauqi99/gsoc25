cmake_minimum_required(VERSION 3.10)
project(GeomScaleGSoC25)

set(CMAKE_CXX_STANDARD 14)

# Locate dependencies
find_package(Eigen3 REQUIRED)
find_package(Volesti REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR} ${Volesti_INCLUDE_DIRS})

