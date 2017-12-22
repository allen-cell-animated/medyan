# Project files for Visual Studio 2017

## Usage

### Main project
+ Configure BOOST include absolute path in `AdditionalIncludeDirectories` in MEDYAN.vcxproj.
+ For debugging in VS, you need to configure the command line arguments in MEDYAN.vcxproj.user.
+ Load MEDYAN.sln to build and debug the MEDYAN project.

### Test project
+ Configure BOOST include absolute path in `AdditionalIncludeDirectories` in MEDYAN_TEST.vcxproj.
+ Configure googletest include absolute path in `AdditionalIncludeDirectories` in MEDYAN_TEST.vcxproj.
+ Configure googletest library absolute path in `AdditionalLibraryDirectories` in MEDYAN_TEST.vcxproj.
+ Change the name of gtest.lib and gtest_main.lib in `AdditionalDependencies` in MEDYAN_TEST.vcxproj if you have different names for those libraries.
+ Load MEDYAN.sln to build MEDYAN_TEST project.

## When new directories are added:
+ Make sure the directories are added to `AdditionalIncludeDirectories` in both MEDYAN.vcxproj and MEDYAN_TEST.vcxproj if you want to avoid paths in the `#include` preprocessor.
