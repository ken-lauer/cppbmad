#!/bin/bash

if [ -z "$ACC_ROOT_DIR" ]; then
  echo "ACC_ROOT_DIR is unset."
  exit 1
fi

build_type=${1:-debug}

python -m codegen &&
  cmake -DCMAKE_BUILD_TYPE="${build_type}" -B "${build_type}" . &&
  make -j 8 -C "${build_type}" &&
  ./$build_type/test_all_encompassing &&
  ./$build_type/test_integration &&
  ./$build_type/test_arrays &&
  cd python/examples/ && python csr.py
