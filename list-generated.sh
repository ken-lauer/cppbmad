#!/bin/bash
# List all generated files, one per line, from the manifest emitted by
# `python -m codegen`. Used by commit-generated.sh, revert-generated.sh,
# .check_untracked.sh, and CI workflows.
grep -v -e '^#' -e '^[[:space:]]*$' generated-manifest.txt
