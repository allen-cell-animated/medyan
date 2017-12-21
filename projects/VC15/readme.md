# Project files for Visual Studio 2017

## Usage:
+ Configure BOOST library absolute path in `AdditionalIncludeDirectories` in MEDYAN/MEDYAN.vcxproj
+ For debugging in VS, you need to configure the command line arguments in project MEDYAN's property configuration.
+ Load MEDYAN.sln to compile and debug.

## When new directories are added:
+ Make sure the directories are added to `AdditionalIncludeDirectories` if you want to avoid paths in the `#include` preprocessor.
