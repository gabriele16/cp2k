#!/bin/bash -e

# Download the repository pin file
# Set noninteractive mode for apt
! export DEBIAN_FRONTEND=noninteractive

# Download the repository pin file
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin
sudo mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600

# Add NVIDIA’s GPG key (the warning about apt-key is okay)
sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub

# Add the CUDA repository non-interactively (the -y flag skips the prompt)
sudo add-apt-repository -y "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/ /"

# Update package lists
sudo apt-get update

# Install CUDA 11.8
sudo apt-get install -y cuda-11-8

# Update environment variables for this session
export PATH=/usr/local/cuda-11.8/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-11.8/lib64:$LD_LIBRARY_PATH



sudo rm -rf /usr/local/cuda
sudo rm -rf /usr/local/cuda-12.5/
sudo rm -rf /usr/local/cuda-12
sudo ln -s /usr/local/cuda-11.8 /usr/local/cuda
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
nvcc --version
