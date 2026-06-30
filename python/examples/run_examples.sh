#!/bin/bash

set -xe

export MPLBACKEND=Agg # no plot windows, please

python bbu.py
python wall_generator.py "$ACC_ROOT_DIR"/bsim/wall_generator/test/lat.bmad 8 0.01
python lat_ele_locator.py
python csr.py

python -c 'import matplotlib, pytao' && {
  python tao-plot.py
  python tao-subproc.py
  python floor.py
} || echo "No matplotlib/PyTao. Skipping those examples."
