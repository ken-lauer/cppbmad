#!/bin/bash

if [ -z "$ACC_ROOT_DIR" ]; then
  echo "ACC_ROOT_DIR is unset."
  exit 1
fi

build_type=${1:-debug}

python -m codegen &&
  cmake -DCMAKE_CXX_FLAGS="-ftime-trace" -DCMAKE_BUILD_TYPE="${build_type}" -B "${build_type}" . &&
  make -j 8 -C "${build_type}" &&
  ./build/test_all_encompassing &&
  ./build/test_integration &&
  ./build/test_arrays &&
  cd python/examples/ && python csr.py
