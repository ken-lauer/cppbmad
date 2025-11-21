#!/bin/bash

if [ -z "$ACC_ROOT_DIR" ]; then
  echo "ACC_ROOT_DIR is unset."
  exit 1
fi

python3 -m codegen &&
  cmake -B build . &&
  make -C build
