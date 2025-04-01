#!/bin/bash -e

#git clone --recursive -b colab-nequip-2024.2 https://github.com/gabriele16/cp2k.git cp2k

# Default option: do not compile CP2K from source.
INSTALL_CP2K="no"

# Parse arguments
for arg in "$@"; do
    case $arg in
        --install-cp2k=*)
            INSTALL_CP2K="${arg#*=}"
            shift
            ;;
        *)
            # unknown option
            ;;
    esac
done


# Print start message
echo "Starting installation..."

# Force reinstall dependencies and log output
echo "Installing nequip and torch..."

pip install wandb
pip install mkl mkl-include
pip install --force-reinstall nequip==0.6 torch==1.13
pip install --force-reinstall numpy==1.26.4
pip install jupyter jupyter_contrib_nbextensions nglview
#jupyter-nbextension enable nglview --py --sys-prefix

# Skip cloning allegro if it already exists
if [ ! -d "/content/allegro" ]; then
  echo "Cloning allegro repository..."
  cd /content && git clone --depth 1 https://github.com/mir-group/allegro.git
else
  echo "allegro directory already exists, skipping clone."
fi

# Force reinstall allegro
echo "Reinstalling allegro..."
pip install --force-reinstall allegro

# Skip downloading libtorch if it already exists
if [ ! -d "/content/libtorch" ]; then
  echo "Downloading libtorch..."
  cd /content && wget https://download.pytorch.org/libtorch/cu118/libtorch-cxx11-abi-shared-with-deps-2.0.0%2Bcu118.zip -O libtorch.zip
  echo "Unzipping libtorch..."
  unzip /content/libtorch.zip -d /content/
else
  echo "libtorch directory already exists, skipping download."
fi
cd /content/libtorch/lib && ln -s libnvrtc-builtins-7237cb5d.so.11.7 libnvrtc-builtins.so.11.8

if [[ "$INSTALL_CP2K" == "yes" ]]; then
    echo "Compiling CP2K from source..."
    # Make sure the environment variables are set so that CUDA is found.
    export PATH=/usr/local/cuda/bin${PATH:+:${PATH}}
    export LD_LIBRARY_PATH=/usr/local/cuda/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

    # Install CP2K toolchain and compile CP2K
    cd cp2k/tools/toolchain
    ./install_cp2k_toolchain.sh \
      --with-spglib=no \
      --with-cosma=no \
      --with-libxsmm=no \
      --with-libvori=no \
      --with-libxc=no \
      --with-plumed=no \
      --with-sirius=no \
      --with-elpa=no \
      --with-libint=no \
      --with-fftw=no \
      --with-libvdwxc=no \
      --with-gsl=no \
      --with-openmpi=install \
      --with-libtorch=/content/libtorch \
      --with-dftd4=no \
      --with-libgrpp=no \
      --with-hdf5=no \
      --with-cusolvermp=no \
      --with-spla=no \
      --with-ninja=no \
      --enable-cuda \
      --gpu-ver=P100

    # Copy architecture files and build CP2K
    cp /content/cp2k/tools/toolchain/install/arch/* /content/cp2k/arch/
    source /content/cp2k/tools/toolchain/install/setup
    cd /content/cp2k/
    make -j 2 ARCH=local_cuda VERSION="ssmp"

    # Prepare prebuilt directory to archive
    echo "Creating prebuilt package..."
    cd /content
    mkdir -p cp2k_prebuilt_cuda/exe
    cp -r cp2k/exe/local_cuda cp2k_prebuilt_cuda/exe/.
    cp -r cp2k/tools/toolchain cp2k_prebuilt_cuda/.
    cp -r cp2k/arch cp2k_prebuilt_cuda/.

    # Create compressed archive of the built CP2K package
    tar -czvf cp2k_prebuilt_cuda.tar.gz cp2k_prebuilt_cuda

    # (Optional) Copy the archive to Google Drive if your drive is mounted.
    if [ -d "/content/drive/MyDrive" ]; then
        cp cp2k_prebuilt_cuda.tar.gz /content/drive/MyDrive/
        echo "Archive copied to Google Drive (/content/drive/MyDrive/)"
    else
        echo "Google Drive not mounted. Archive is in the current directory."
    fi

else
    echo "Downloading precompiled CP2K package..."
    # Install gdown if not already installed
    cd /content
    pip install gdown

    # Download the precompiled CP2K package from Google Drive using its file ID.
    # Replace the file ID below with the correct one for your precompiled package.
    
    gdown https://drive.google.com/uc?id=1fVG4Fm92XN_vJNGcgkabu7unLM5sqbRm

    # Unpack the precompiled package
    tar -xzvf cp2k_prebuilt_cuda.tar.gz

    # Copy the executable and architecture files into the CP2K source directory.
    cp -r cp2k_prebuilt_cuda/exe /content/cp2k/.
    cp -r cp2k_prebuilt_cuda/arch/* /content/cp2k/arch/.

    # Rename the existing toolchain directory and replace it with the precompiled one.
    mv /content/cp2k/tools/toolchain /content/cp2k/tools/toolchain_not_built
    cp -r cp2k_prebuilt_cuda/toolchain /content/cp2k/tools/.
fi

