#!/bin/sh -e

# The build script for Deepthought2 HPC cluster

medyan_root_dir=$(X= cd -- "$(dirname -- "$0")" && pwd -P)

# Set up variables
export MEDYAN_BOOST_INSTALL_MODE="manual"
export MEDYAN_BOOST_INCLUDE_DIR="$(BOOST_INCLUDE)"
export MEDYAN_BOOST_LIBRARY_DIR="$(BOOST_LIB)"
export MEDYAN_BOOST_LIBRARIES="boost_system"
export MEDYAN_ADDITIONAL_LINK_DIRS="$(GCC_ROOTDIR)/lib64"
export MEDYAN_RPATH="$(GCC_ROOTDIR)/lib64"

# Run the script
$medyan_root_dir/build-medyan.sh
