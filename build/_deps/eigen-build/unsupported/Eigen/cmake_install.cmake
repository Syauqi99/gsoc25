# Install script for directory: /home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/AdolcForward"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/AlignedVector3"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/ArpackSupport"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/AutoDiff"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/BVH"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/EulerAngles"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/FFT"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/IterativeSolvers"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/KroneckerProduct"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/LevenbergMarquardt"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/MatrixFunctions"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/MoreVectorization"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/MPRealSupport"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/NonLinearOptimization"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/NumericalDiff"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/OpenGLSupport"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/Polynomials"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/Skyline"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/SparseExtra"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/SpecialFunctions"
    "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/Splines"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/home/syauqirp/gsoc25/build/_deps/eigen-src/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/syauqirp/gsoc25/build/_deps/eigen-build/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

