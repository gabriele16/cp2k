#!/bin/bash
# Usage: ./run_cp2k.sh -i input_file

# Parse command line options
while getopts "i:" opt; do
  case "$opt" in
    i) input_file="$OPTARG" ;;
    *) echo "Usage: $0 -i input_file" >&2; exit 1 ;;
  esac
done

# Check that an input file was provided
if [ -z "$input_file" ]; then
  echo "Usage: $0 -i input_file" >&2
  exit 1
fi

# Set LD_LIBRARY_PATH so that CP2K finds the proper libraries
export LD_LIBRARY_PATH=/content/libtorch/lib:/usr/local/cuda/lib64:$LD_LIBRARY_PATH

# Run CP2K with the provided input file
/content/cp2k/exe/local_cuda/cp2k.ssmp -i "$input_file"

