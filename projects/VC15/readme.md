# Project files for Visual Studio 2017

## Usage

+ Configure BOOST include absolute path in `AdditionalIncludeDirectories` in MEDYAN.vcxproj.

### Main project
+ For debugging in VS, you need to configure the command line arguments in MEDYAN.vcxproj.user.
+ Load MEDYAN.sln, choose `Debug` or `Release` configuration to build the MEDYAN project.

### Test project
+ Configure googletest include absolute path in `AdditionalIncludeDirectories` in MEDYAN.vcxproj (in `TestDebug` or `TestRelease` configurations only).
+ Configure googletest library absolute path in `AdditionalLibraryDirectories` in MEDYAN.vcxproj (in `TestDebug` or `TestRelease` configurations only).
+ Change the name of gtest.lib and gtest_main.lib in `AdditionalDependencies` in MEDYAN.vcxproj if you have different names for those libraries.
+ Load MEDYAN.sln, choose `TestDebug` or `TestRelease` configuration to build MEDYAN project.

## When new directories are added:
+ Make sure the directories are added to `AdditionalIncludeDirectories` in MEDYAN.vcxproj if you want to avoid relative paths in the `#include` preprocessor.
