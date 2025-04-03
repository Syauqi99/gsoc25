#!/bin/bash

# Create third_party directory if it doesn't exist
mkdir -p third_party

# Download Eigen if it doesn't exist
if [ ! -d "third_party/eigen" ]; then
    echo "Downloading Eigen..."
    cd third_party
    git clone --depth 1 https://gitlab.com/libeigen/eigen.git
    cd ..
    echo "Eigen downloaded successfully."
else
    echo "Eigen already exists in third_party/eigen."
fi

# Make script executable
chmod +x download_eigen.sh
