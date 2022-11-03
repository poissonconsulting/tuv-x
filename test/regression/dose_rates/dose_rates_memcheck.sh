#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# copy photolysis rate configuration into template
sed -e '/#DOSE_RATES/ {' -e 'r data/dose_rates.json' -e 'd' -e '}' -i test/data/doserates.test.config.json

exec_oldtuv() {
  ./oldtuv DO_RAYLEIGH DO_O2 DO_O3 DO_AEROSOLS DO_CLOUDS < test/regression/tuv_scenario_2.in
}
exec_newtuv() {
  valgrind --error-exitcode=1 --trace-children=yes --leak-check=full --gen-suppressions=all --suppressions=test/valgrind.supp ./tuv-x test/data/doserates.test.config.json
}
exec_analysis() {
  python3 test/regression/dose_rates/sw.compare.py test/regression/dose_rates odat/OUTPUTS output
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
