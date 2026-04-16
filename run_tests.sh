#!/bin/bash

build_type=${1:-debug}

./$build_type/cppbmad_tests &&
  (cd python/examples/ && bash run_examples.sh) &&
  (cd python/tests/ && python -m pytest)
