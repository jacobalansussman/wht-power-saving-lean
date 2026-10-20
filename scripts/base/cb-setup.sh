#!/bin/bash
# c-base step 1-2a: clone Shea's repository at the proof commit, fetch Mathlib d13f23b and its cache.
. /content/work/env.sh
set -x
date -u
cd /content/work
git clone https://github.com/shea256/fourier-transform-below-nlogn /content/work/shea || exit 11
cd /content/work/shea
git checkout --detach 5239484 || exit 12
git rev-parse HEAD
git log -1 --format='%H %ci %s'
git rev-parse origin/HEAD origin/main 2>/dev/null
git ls-remote origin | head -5
git diff --stat 5239484 0e61af6 -- lean/ | tail -3
echo "diff-lean-lines: $(git diff 5239484 0e61af6 -- lean/ | wc -l)"
cd lean/complex588
cat lean-toolchain
grep -n '"rev"' lake-manifest.json | head -2
/usr/bin/time -v -o /content/work/logs/cb-cache.time lake exe cache get
echo "CACHE rc=$?"
tail -20 /content/work/logs/cb-cache.time
lean --version; lake --version
du -sh .lake/packages $MATHLIB_CACHE_DIR /content/work/elan
date -u
