#!/bin/bash -e

# Print start message
echo "Starting installation..."

# Force reinstall dependencies and log output
echo "Installing nequip and torch..."
pip install --force-reinstall nequip=0.6 torch==1.13

# Skip cloning cp2k if it already exists
if [ ! -d "/content/cp2k" ]; then
  echo "Cloning cp2k repository..."
  cd /content && git clone --depth 1 https://github.com/cp2k/cp2k.git
else
  echo "cp2k directory already exists, skipping clone."
fi

# Skip cloning allegro if it already exists
if [ ! -d "/content/allegro" ]; then
  echo "Cloning allegro repository..."
  cd /content && git clone --depth 1 https://github.com/mir-group/allegro.git
else
  echo "allegro directory already exists, skipping clone."
fi

# Force reinstall allegro
echo "Reinstalling allegro..."
pip install --force-reinstall allegro/

# Skip downloading libtorch if it already exists
if [ ! -d "/content/libtorch" ]; then
  echo "Downloading libtorch..."
  cd /content && wget https://download.pytorch.org/libtorch/cu118/libtorch-cxx11-abi-shared-with-deps-2.0.0%2Bcu118.zip -O libtorch.zip
  echo "Unzipping libtorch..."
  unzip /content/libtorch.zip -d /content/
else
  echo "libtorch directory already exists, skipping download."
fi

