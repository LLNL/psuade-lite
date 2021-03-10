#!/bin/bash

set -x -e

cmake "$SRC_DIR" -DCMAKE_CXX_FLAGS="--std=c++14" -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_OSX_SYSROOT="$CONDA_BUILD_SYSROOT"
cmake --build . --target install
