#!/usr/bin/env bash
# run-tests.sh — execute the JWST-L2 pytest suite.
# Usage:
#   ./run-tests.sh           # run all tests verbose
#   ./run-tests.sh -k geom   # filter to tests matching keyword
set -e
cd "$(dirname "$0")"
# --basetemp pinned to avoid a pytest cleanup-recursion edge case on bind-mounted
# filesystems (observed in Cowork sandbox; harmless on native macOS/Linux but
# pinning makes the behaviour identical across environments).
python3 -m pytest tests/ -v --tb=short -p no:cacheprovider --basetemp=/tmp/pytest-jwst-l2 "$@"
