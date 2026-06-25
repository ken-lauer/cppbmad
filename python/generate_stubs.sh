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

# Post-processing workaround for https://github.com/sizmailov/pybind11-stubgen/pull/272 :
strip_self_prefix() {
  local f=$1
  if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' "s/${MODULE_NAME}\.//g" "$f"
  else
    sed -i "s/${MODULE_NAME}\.//g" "$f"
  fi
}

if [[ -f "$OUTPUT_DIR/$MODULE_NAME/__init__.pyi" ]]; then
  strip_self_prefix "$OUTPUT_DIR/$MODULE_NAME/__init__.pyi"
else
  echo "ERROR: no stub output found for $MODULE_NAME under $OUTPUT_DIR" >&2
  exit 1
fi

echo "Stub generation complete."
