#!/bin/sh -e

# Set directories and files
medyan_root_dir=$(X= cd -- "$(dirname -- "$0")" && pwd -P)
medyan_build_dir="$medyan_root_dir/build"

medyan_vcpkg_cmake_toolchain="$medyan_root_dir/scripts/.build/vcpkg/scripts/buildsystems/vcpkg.cmake"

# Set variables
if [ "$MEDYAN_BOOST_INSTALL_MODE" = "manual" ]; then
    medyan_cmake_boost="-DMEDYAN_BOOST_INSTALL_MODE=manual"
    medyan_cmake_boost="$medyan_cmake_boost -DMEDYAN_BOOST_INCLUDE_DIR=$MEDYAN_BOOST_INCLUDE_DIR"
    medyan_cmake_boost="$medyan_cmake_boost -DMEDYAN_BOOST_LIBRARY_DIR=$MEDYAN_BOOST_LIBRARY_DIR"
    medyan_cmake_boost="$medyan_cmake_boost -DMEDYAN_BOOST_LIBRARIES=$MEDYAN_BOOST_LIBRARIES"
    medyan_need_install_boost=false
elif [ "$MEDYAN_BOOST_INSTALL_MODE" = "find" ]; then
    medyan_cmake_boost="-DMEDYAN_BOOST_INSTALL_MODE=find"
    medyan_need_install_boost=false
else
    medyan_need_install_boost=true
fi
if [ -n "$MEDYAN_ADDITIONAL_LINK_DIRS" ]; then
    medyan_cmake_additional_link_dirs="-DMEDYAN_ADDITIONAL_LINK_DIRS=$MEDYAN_ADDITIONAL_LINK_DIRS"
fi
if [ -n "$MEDYAN_RPATH" ]; then
    medyan_cmake_rpath="-DMEDYAN_RPATH=$MEDYAN_RPATH"
fi

export medyan_root_dir
export medyan_need_install_boost

# Run bootstrapping
echo "Bootstrapping..."
$medyan_root_dir/scripts/bootstrap.sh

# Create build files
echo "Creating build files..."
mkdir -p $medyan_build_dir
cd $medyan_build_dir
cmake \
    $medyan_cmake_boost \
    $medyan_cmake_additional_link_dirs \
    $medyan_cmake_rpath \
    .. $medyan_vcpkg_cmake_toolchain

# Build medyan
echo "Building medyan..."
if [ -n "$MEDYAN_BUILD_PROCS" ]; then
    medyan_make_j="-j$MEDYAN_BUILD_PROCS"
fi
make $medyan_make_j
