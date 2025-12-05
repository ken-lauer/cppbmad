#!/bin/bash

if [ -z "$ACC_ROOT_DIR" ]; then
  echo "ACC_ROOT_DIR is unset."
  exit 1
fi

python -m codegen &&
  cmake -B build . &&
  make -j -C build &&
  ./build/test_all_encompassing &&
  ./build/test_integration
