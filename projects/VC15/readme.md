# Project files for Visual Studio 2017

## Usage:
+ Configure BOOST library absolute path in `AdditionalIncludeDirectories` in MEDYAN/MEDYAN.vcxproj
+ Load MEDYAN.sln to compile and debug.

## When new directories are added:
+ Make sure the directories are added to `AdditionalIncludeDirectories` if you want to avoid paths in the `#include` preprocessor.
