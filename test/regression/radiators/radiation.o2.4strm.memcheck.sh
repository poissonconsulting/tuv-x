#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v

exec_oldtuv() {
  ./oldtuv DO_O2 < test/regression/tuv_scenario_6.in
}
exec_newtuv() {
  valgrind --error-exitcode=1 --trace-children=yes --leak-check=full --gen-suppressions=all --suppressions=test/valgrind.supp ./tuv-x test/data/radiators.o2.4strm.config.json
}
exec_analysis() {
  python3 tool/diagnostics/var.compare.py test/regression/radiators/radiation.o2.compare.json
}

if ! exec_oldtuv; then
  echo FAIL - old TUV
  exit 1
fi

if ! exec_newtuv; then
  echo FAIL - new TUV
  exit 1
fi

if ! exec_analysis; then
  echo FAIL - analysis
  exit 1
fi

echo PASS
exit 0
