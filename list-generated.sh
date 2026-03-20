#!/bin/bash
# List all generated files tracked by git.
# Used by commit-generated.sh, revert-generated.sh, and CI workflows.
echo $(git ls-files | grep generated) \
  python/pybmad/_*.pyi \
  python/pybmad/__init__.py \
  python/pybmad/_*.py \
  coverage.html \
  python/docs/source/api/*.md
