#!/bin/bash
# fs-build: setup -> build+checks -> official comparator, each only if the previous one ended well.
L=/content/work/logs/fs-build; F=/content/work/fs-build
until grep -q "^EXIT" $L/setup.out 2>/dev/null; do sleep 15; done
grep -q "^EXIT 0" $L/setup.out && [ -e $F/SETUP-OK ] || { echo "setup failed: $(tail -3 $L/setup.out)"; exit 1; }
grep -q "lake --no-build U8XM rc=0 : All targets up-to-date" $L/setup.out || { echo "U8XM not up to date in my copy: stop"; exit 2; }
grep -q "all identical" $L/setup.out && grep -q "U8XMHH-SOURCES-OK" $L/setup.out && grep -q "XHOST-SOURCES-OK" $L/setup.out || { echo "hash checks missing: stop"; exit 3; }
bash $F/scripts/fsb-build.sh > $L/build.out 2>&1; rc=$?; echo "build rc=$rc $(date -u)"
[ $rc -eq 0 ] && [ -e $F/BUILD-CHECKS-DONE ] || exit 4
grep -q "RESULT: PASS" $L/compare-pos.log && grep -q "RESULT: FAIL" $L/compare-neg.log || { echo "Compare.lean results not as required: no comparator"; exit 5; }
bash $F/scripts/fsb-cmp.sh > $L/cmp.out 2>&1; echo "cmp rc=$? $(date -u)"
