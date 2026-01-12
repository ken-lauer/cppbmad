#!/usr/bin/env bash

paths_to_check=(
  "src/generated/*.cpp"
  "include/bmad/generated/*.hpp"
  "python/src/generated/*.cpp"
  "python/include/pybmad/generated/*.hpp"
)

untracked_files=""
for path in "${paths_to_check[@]}"; do
  found=$(git ls-files --others --exclude-standard "$path")
  if [ -n "$found" ]; then
    untracked_files="${untracked_files}${found}\n"
  fi
done

if [ -n "$untracked_files" ]; then
  echo "Error: Untracked files found:" >/dev/stderr
  printf "%b" "$untracked_files"
  echo "To add these:" >/dev/stderr
  echo "git add \$(bash .check_untracked.sh)" >/dev/stderr

  exit 1
fi
