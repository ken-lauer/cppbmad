#!/usr/bin/env bash
# Fail if any file in the generated manifest is untracked by git.
# Run after `python -m codegen` to catch newly generated files that
# nobody remembered to `git add`.

untracked_files=""
while IFS= read -r path; do
  [ -z "$path" ] && continue
  found=$(git ls-files --others --exclude-standard "$path")
  if [ -n "$found" ]; then
    untracked_files="${untracked_files}${found}\n"
  fi
done < <(bash list-generated.sh)

if [ -n "$untracked_files" ]; then
  echo "Error: Untracked generated files found:" >/dev/stderr
  printf "%b" "$untracked_files"
  echo "To add these:" >/dev/stderr
  echo "git add \$(bash list-generated.sh)" >/dev/stderr

  exit 1
fi
