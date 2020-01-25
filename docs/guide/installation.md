# Installation

To use MEDYAN, one needs to build it from the source code. The source code is available for download at <http://medyan.org/download.html>. One can also clone the git repository for the latest version.

## Requirements

A C++14 compatible compiler such as the following is required:

- gcc 5 and above
- clang 3.4 and above
- MSVC with Visual Studio 2017 15.0 and above

Additionally, the following external library is required:

- Boost libraries 1.49 and above

The Boost library can be downloaded at <https://www.boost.org/users/download/>.

## Building steps

MEDYAN provides multiple build configuration files and can be compiled on most of the platforms. For example, the `Makefile` is provided to be used with `make`. On MacOS, an Xcode project file is provided. On Windows, a Visual Studio solution is also available.

### Build with `make`

This section introduces building MEDYAN using `make`. The `Makefile` is available in `src` directory. Before running `make`, several changes are necessary in `src/Makefile`.

- If `g++` is not in your path, then replace `g++` in `CXX = g++ ...` with the path to your actual C++ compiler.
- In `INCLUDES = ...`, add the path to the boost library as `-I<your-boost-install-location>`.

A number of predefined macros are also provided in the `Makefile`. In general, you are not expected to change any one of them. To view all the predefined macros and their usage, refer to [Predefined macros](../manual/predefined-macro.md).

Finally, make sure you are in `src` directory and run the following command:
```console
$ make
```
to start building. If all the configurations are correct, an executable named `MEDYAN` should be created in the current directory.

In the `make` command, one can use `-j` option as `make -j<number-of-tasks>` to utilize multiple cores to compile parallelly for faster compiling.

### Build with Xcode (MacOS)

An Xcode project file is available in the root directory of the source code. On MacOS, it can be used with Xcode to build, run and debug the MEDYAN program.

To build the program, simply load `MEDYAN.xcodeproj` in Xcode.

**TODO configurations**

Finally, click the "build and run" button to start building.

### Build with Visual Studio (Windows)

A Visual Studio solution file is provide in `projects/VS16` directory, and can be used with Visual Studio to build, run and debug the MEDYAN program.

To build, simply load `MEDYAN.sln` into Visual Studio.

**TODO configuration**

Finally, one can start to build and run the program.
