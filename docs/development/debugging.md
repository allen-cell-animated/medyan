# Debugging MEDYAN

## VS Code

# Build MEDYAN for debug

```console
> MEDYAN_CMAKE_EXTRA_ARGS="-DCMAKE_BUILD_TYPE=Debug" ./conf.sh
> cd build
> make
```

# Set command line arguments

1. Open `.vscode/launch.json`.

2. edit the `"args": ["test"],` line to include the command line arguments you want to debug with.

2. edit the `"cwd": "${workspaceFolder}/build",` line to choose the directory to run MEDYAN from.