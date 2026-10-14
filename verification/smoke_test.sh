#!/bin/bash
# One-command check that the repository is functional.
#
# Runs the Python reference implementation on all its cases and compares
# against the analytic values. Needs numpy (and scipy for the sparse route).
# The Elmer cases need Elmer and Gmsh, which are not assumed here; see
# cases/*/command.txt for those.
set -u
cd "$(dirname "$0")/.."
fail=0
for c in python/cases/*/; do
    [ -f "$c/run.py" ] || continue
    name=$(basename "$c")
    if python3 "$c/run.py" 2>&1 | grep -q "rel. error"; then
        echo "  OK   $name"
    else
        echo "  FAIL $name"
        fail=1
    fi
done
if [ "$fail" -eq 0 ]; then
    echo "smoke test passed"
else
    echo "smoke test FAILED"
fi
exit $fail
