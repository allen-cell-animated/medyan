# Project files for Visual Studio 2017

## Usage

### Main project
+ Configure BOOST library absolute path in `AdditionalIncludeDirectories` in MEDYAN.vcxproj
+ For debugging in VS, you need to configure the command line arguments in MEDYAN.vcxproj.user.
+ Load MEDYAN.sln to build and debug the MEDYAN project.

## When new directories are added:
+ Make sure the directories are added to `AdditionalIncludeDirectories` if you want to avoid paths in the `#include` preprocessor.
