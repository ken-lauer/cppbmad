#!/usr/bin/env bash

set -e

MODULE_NAME=${1-_pybmad}
OUTPUT_DIR=${2-./pybmad}
STUBGEN_EXECUTABLE=${3-pybind11-stubgen}

if ! python -c "import _pybmad" 2>/dev/null; then
  export PYTHONPATH=$PWD/pybmad:$PYTHONPATH
fi

echo "Generating stubs for: $MODULE_NAME"
echo "Output directory: $OUTPUT_DIR"

"$STUBGEN_EXECUTABLE" -o "$OUTPUT_DIR" "$MODULE_NAME" --ignore-all-errors

# Post-processing: Replace "_pybmad." with ""
STUB_FILE=$OUTPUT_DIR/$MODULE_NAME.pyi

# actual fix pending: https://github.com/sizmailov/pybind11-stubgen/pull/272

if [[ "$OSTYPE" == "darwin"* ]]; then
  sed -i '' 's/_pybmad\.//g' "$STUB_FILE"
else
  sed -i 's/_pybmad\.//g' "$STUB_FILE"
fi

echo "Stub generation complete."
