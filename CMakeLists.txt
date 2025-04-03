# GSoC 2025 Project
# Based on VolEsti (volume computation and sampling library)

project(GSoC2025)

CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

# Locate Intel MKL root (in case it is enabled)
if (APPLE)
    set(MKLROOT /opt/intel/oneapi/mkl/latest)
elseif(UNIX)
    set(MKLROOT /opt/intel/oneapi/mkl/latest)
endif()

# Options
option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)
option(BUILTIN_EIGEN "Use eigen from volesti/external" OFF)
option(BUILTIN_AUTODIFF "Use autodiff from volesti/external" ON)
option(USE_MKL "Use MKL library to build eigen" OFF)

if(DISABLE_NLP_ORACLES)
  add_definitions(-DDISABLE_NLP_ORACLES)
else()
  find_library(IFOPT NAMES libifopt_core.so PATHS /usr/local/lib)
  find_library(IFOPT_IPOPT NAMES libifopt_ipopt.so PATHS /usr/local/lib)
  find_library(GMP NAMES libgmp.so PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)
  find_library(MPSOLVE NAMES libmps.so PATHS /usr/local/lib)
  find_library(PTHREAD NAMES libpthread.so PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)
  find_library(FFTW3 NAMES libfftw3.so.3 PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)

  if (NOT IFOPT)
    message(FATAL_ERROR "This program requires the ifopt library, and will not be compiled.")
  elseif (NOT GMP)
    message(FATAL_ERROR "This program requires the gmp library, and will not be compiled.")
  elseif (NOT MPSOLVE)
    message(FATAL_ERROR "This program requires the mpsolve library, and will not be compiled.")
  elseif (NOT FFTW3)
    message(FATAL_ERROR "This program requires the fftw3 library, and will not be compiled.")
  else()
    message(STATUS "Library ifopt found: ${IFOPT}")
    message(STATUS "Library gmp found: ${GMP}")
    message(STATUS "Library mpsolve found: ${MPSOLVE}")
    message(STATUS "Library fftw3 found: ${FFTW3}")
  endif()
endif()

# Set CMake policy
if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Include path to volesti's cmake files
set(VOLESTI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tp_library/volesti)

# Include the dependencies using volesti's cmake files
include("${VOLESTI_DIR}/external/cmake-files/Autodiff.cmake")
GetAutodiff()

include("${VOLESTI_DIR}/external/cmake-files/Eigen.cmake")
GetEigen()

include("${VOLESTI_DIR}/external/cmake-files/Boost.cmake")
GetBoost()

include("${VOLESTI_DIR}/external/cmake-files/LPSolve.cmake")
GetLPSolve()

include("${VOLESTI_DIR}/external/cmake-files/QD.cmake")
GetQD()

# Code Coverage Configuration
add_library(coverage_config INTERFACE)

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O1        # O0 (or no) optimization takes too much time and causes CircleCI test failure.
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif()

if (BUILTIN_AUTODIFF)
    include_directories(BEFORE ${VOLESTI_DIR}/external/_deps/Autodiff)
else()
    include_directories(BEFORE /usr/local/include)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

if (USE_MKL)
  find_library(BLAS NAMES libblas.so libblas.dylib PATHS /usr/local/Cellar/lapack/3.9.1_1/lib /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu /usr/local/Cellar/openblas/0.3.15_1/lib /usr/lib)
  find_library(GFORTRAN NAME libgfortran.dylib PATHS /usr/local/Cellar/gcc/10.2.0_4/lib/gcc/10)
  find_library(LAPACK NAME liblapack.dylib PATHS /usr/lib)
  find_library(OPENMP NAME libiomp5.dylib PATHS /opt/intel/oneapi/compiler/2021.1.1/mac/compiler/lib)
  include_directories(BEFORE ${MKLROOT}/include)
  set(PROJECT_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GFORTRAN_LIBRARIES})
  set(MKL_LINK "-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
  add_definitions(-DEIGEN_USE_MKL_ALL)
else()
  set(MKL_LINK "")
endif()

# Include volesti directories
include_directories(BEFORE ${VOLESTI_DIR}/external)
include_directories(BEFORE ${VOLESTI_DIR}/include)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# For Eigen
if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
  add_compile_options(-D "EIGEN_NO_DEBUG")
else()
  add_compile_definitions("EIGEN_NO_DEBUG")
endif()

# Compiler flags
add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
set(CMAKE_CXX_STANDARD 17)  # enable the c++17 support needed by autodiff
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")
add_definitions(${CMAKE_CXX_FLAGS} "-DMKL_ILP64")

# Add your executable
add_executable(my_program scripts/test2_main.cpp)

# Link against volesti dependencies
TARGET_LINK_LIBRARIES(my_program lp_solve ${MKL_LINK} coverage_config)

# If you need NLP oracles
if(NOT DISABLE_NLP_ORACLES)
  TARGET_LINK_LIBRARIES(my_program ${IFOPT} ${IFOPT_IPOPT} ${PTHREAD} ${GMP} ${MPSOLVE} ${FFTW3})
endif()

# Add additional executables as needed
# add_executable(another_program src/another.cpp)
# TARGET_LINK_LIBRARIES(another_program lp_solve ${MKL_LINK} coverage_config)

# Enable testing if needed
# enable_testing()
# add_test(NAME test_example COMMAND my_program --test)

# Set additional flags for specific executables if needed
set(ADDITIONAL_FLAGS "-march=native -DSIMD_LEN=0 -DTIME_KEEPING")
# set_target_properties(my_program PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS})

# Add new project LinearProgramming
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(InteriorPointLP VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add Eigen dependency
find_package(Eigen3 REQUIRED)

# Include directories
include_directories(${EIGEN3_INCLUDE_DIR})

# Add the source files
set(SOURCES
    src/interior_point_lp.cpp
    src/lp_utils.cpp
)

# Add the header files
set(HEADERS
    include/interior_point_lp.h
    include/lp_utils.h
)

# Create the library
add_library(interior_point_lp STATIC ${SOURCES} ${HEADERS})

# Link Eigen to the library
target_link_libraries(interior_point_lp Eigen3::Eigen)

# Add executable for testing
add_executable(test3_main scripts/test3_main.cpp)

# Link the interior_point_lp library and Eigen to test3_main
target_link_libraries(test3_main interior_point_lp Eigen3::Eigen)

