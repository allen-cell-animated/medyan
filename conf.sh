#!/bin/sh -e

# Set directories and files
export medyan_root_dir=$(X= cd -- "$(dirname -- "$0")" && pwd -P)

# Run bootstrapping
echo "Bootstrapping..."
$medyan_root_dir/scripts/bootstrap.sh
