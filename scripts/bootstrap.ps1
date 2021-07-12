# This script might change the current directory

$ErrorActionPreference = 'Stop'

if(-NOT $medyan_root_dir) {
    Write-Error "Variable medyan_root_dir is not specified"
    exit 1
}

# Set directories
$medyan_build_dir = "$medyan_root_dir\build"
$script_dir = "$medyan_root_dir\scripts"
$build_dir = "$script_dir\.build"
$vcpkg_dir = "$build_dir\vcpkg"

# Set variables
$medyan_vcpkg_cmake_toolchain = "$medyan_root_dir\scripts\.build\vcpkg\scripts\buildsystems\vcpkg.cmake"

if($MEDYAN_ADDITIONAL_LINK_DIRS) {
    $medyan_cmake_additional_link_dirs = "-DMEDYAN_ADDITIONAL_LINK_DIRS=$MEDYAN_ADDITIONAL_LINK_DIRS"
}

if("$MEDYAN_RPATH") {
    $medyan_cmake_rpath = "-DMEDYAN_RPATH=$MEDYAN_RPATH"
}


# Download and setup vcpkg
Function Install-Vcpkg([bool]$required, [bool]$rebuild) {

    if($required) {
        if(-NOT $(Test-Path $vcpkg_dir) -OR $rebuild) {
            Write-Host "Configuring vcpkg..."
            Set-Location $build_dir
            git clone https://github.com/Microsoft/vcpkg.git
            Set-Location $vcpkg_dir
            git checkout 2021.05.12
            .\bootstrap-vcpkg.bat
        } else {
            Write-Host "vcpkg is already installed."
        }
    } else {
        Write-Host "Skipping vcpkg installation."
    }

}

# Setup dependencies
Function Install-VcpkgPackages() {

    $Env:VCPKG_DEFAULT_TRIPLET="x64-windows"

    Set-Location $vcpkg_dir
    .\vcpkg install catch2 eigen3 spectra boost-signals2 boost-pool boost-heap boost-ublas boost-range
}

# Generate using CMake
Function Use-Cmake() {

    Write-Host "Creating build files..."
    mkdir -Force $medyan_build_dir
    Set-Location $medyan_build_dir
    cmake `
        $medyan_cmake_additional_link_dirs `
        $medyan_cmake_rpath `
        .. "-DCMAKE_TOOLCHAIN_FILE=$medyan_vcpkg_cmake_toolchain"
}

# Make directories
mkdir -Force $build_dir

# Use vcpkg to resolve dependencies
Install-Vcpkg $true $false
Install-VcpkgPackages

# Use CMake to generate build files
Use-Cmake
