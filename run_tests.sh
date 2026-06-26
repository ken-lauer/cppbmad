#!/bin/bash

build_type=${1:-debug}

export PYTHONPATH=$PWD/$build_type/python:$PWD/python:$PYTHONPATH

./$build_type/cppbmad_tests &&
  (cd python/examples/ && bash run_examples.sh) &&
  (python -m pytest python/tests/)
