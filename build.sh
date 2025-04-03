#!/bin/bash

# Download Eigen if needed
./download_eigen.sh

# Create build directory if it doesn't exist
mkdir -p build

# Build the project
cd build
cmake ..
make -j4

echo "Build complete. The executable is in build/test3_main"
