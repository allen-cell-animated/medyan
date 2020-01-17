#!/bin/sh

if [ -z $medyan_root_dir ]; then
    echo "Error: variable medyan_root_dir needs to be specified"
    exit 1
fi

# Set directories
cur_dir="$medyan_root_dir/scripts"
build_dir="$cur_dir/.build"
vcpkg_dir="$build_dir/vcpkg"

# Download and setup vcpkg
setup_vcpkg() {
    required=$1
    rebuild=$2
    if [ $required = true ]; then
        if [ ! -d "$vcpkg_dir" -o $rebuild = true ]; then
            echo "Configuring vcpkg..."
            (
                cd $build_dir &&
                git clone https://github.com/Microsoft/vcpkg.git &&
                cd $vcpkg_dir &&
                ./bootstrap-vcpkg.sh
            )
        else
            echo "vcpkg is already installed."
        fi
    else
        echo "Skipping vcpkg installation."
    fi
    return 0
}

# Setup dependencies
vcpkg_install() {
    need_boost=$1

    (
        cd $vcpkg_dir && {
            ./vcpkg install catch2 eigen3 spectra
            if [ "$need_boost" = true ]; then
                ./vcpkg install boost
            fi
        }
    )
}

# Make directories
mkdir -p $build_dir

# Use vcpkg to resolve dependencies
setup_vcpkg true false
vcpkg_install $medyan_need_install_boost
