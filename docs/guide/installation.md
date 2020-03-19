# Installation

To use MEDYAN, one needs to build it from the source code. The source code is available for download at <http://medyan.org/download.html>. One can also clone the git repository for the latest version.

## Requirements

A C++14 compatible compiler such as the following is required:

- gcc >= 5
- clang >= 3.4
- MSVC with Visual Studio 2017 15.0 and above

The following tools are required:

- CMake >= 3.12
- git >= 2.7.0

`medyan` is also dependent on several external libraries, which will be automatically installed before building.

## Building `medyan`

The preferred way of building medyan is to generate the build files using CMake. From the MEDYAN root directory, use either of the following scripts to generate build files.

- `conf.sh` on Linux and MacOS. It will generate `Makefile` in the `build` directory by default.
- `conf.ps1` on Windows. It will generate the Visual Project solution in the `build` directory by default.

Under the hood, the scripts use `vcpkg` to install and configure the required external libraries, and use `cmake` to generate the build files. The external library installation normally happens only for the first build. External variables can be used to control some of the behaviors of the script.

Among all the external libraries required by `medyan`, `boost` is treated a bit differently, simply because `boost` is too large to install quickly, and if the user already has `boost` installed on the computer, it is redundant to do a slow fresh install. One can configure the boost install via the `MEDYAN_BOOST_INSTALL_MODE` variable. The behavior is listed in the following table.

| `MEDYAN_BOOST_INSTALL_MODE` | behavior|
|-----------------------------|---------|
| *not set* (default) | Do a fresh boost install. |
| `manual`            | Manually configure the existing boost library. Additionally, `MEDYAN_BOOST_INCLUDE_DIR`, `MEDYAN_BOOST_LIBRARY_DIR` and `MEDYAN_BOOST_LIBRARIES` are required to locate the boost library. |
| `auto`              | Automatically detect existing boost library on the computer using CMake. The build would fail if CMake cannot find boost. |

On MacOS, one can also change the target build system of the CMake from `make` to Xcode project, by setting `MEDYAN_BUILD_TOOL` to be `"Xcode"`.

**Example 1**: I am building medyan on Linux, and I do not have boost installed anywhere.

```console
> ./conf.sh
> cd build
> make
```

And the medyan would be built in the `build` directory.

**Example 2**: I have boost on my Mac and I want to build medyan.

```console
> MEDYAN_BOOST_INSTALL_MODE="auto" ./conf.sh
```

And an Xcode project file would be in the `build` directory.

**Example 3**: I have boost on my Windows and I do not think CMake can find it anywhere.

```console
> $MEDYAN_BOOST_INSTALL_MODE = "manual"
> $MEDYAN_BOOST_INCLUDE_DIR = "path\to\boost\include"
> $MEDYAN_BOOST_LIBRARY_DIR = "path\to\boost\lib"
> $MEDYAN_BOOST_LIBRARIES = "some_boost_lib_1;some_boost_lib_2"
> .\conf.ps1
```

And a Visual Studio solution file would be in the `build` directory.

Alternatively, medyan provides the `Makefile` to be used with `make`. On MacOS, an Xcode project file is provided. On Windows, a Visual Studio solution is also available. *Note that these are currently maintained for legacy reasons, but the support for these might be dropped in the future.*

### (Old) Build with `make`

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

### (Old) Build with Xcode (MacOS)

An Xcode project file is available in the root directory of the source code. On MacOS, it can be used with Xcode to build, run and debug the MEDYAN program.

To build the program, simply load `MEDYAN.xcodeproj` in Xcode.

Finally, click the "build and run" button to start building.

### (Old) Build with Visual Studio (Windows)

A Visual Studio solution file is provide in `projects/VS16` directory, and can be used with Visual Studio to build, run and debug the MEDYAN program.

To build, simply load `MEDYAN.sln` into Visual Studio.

Finally, one can start to build and run the program.
