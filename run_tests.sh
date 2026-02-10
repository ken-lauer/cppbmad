#!/bin/bash

build_type=${1:-debug}

./$build_type/test_all_encompassing &&
  ./$build_type/test_integration &&
  ./$build_type/test_arrays &&
  (cd python/examples/ && bash run_examples.sh) &&
  (cd python/tests/ && python -m pytest)
