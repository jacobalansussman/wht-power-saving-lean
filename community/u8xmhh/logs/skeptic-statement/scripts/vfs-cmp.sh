#!/bin/bash
# v-fs-statement (skeptic): replay of the official comparator (same binaries as fs-build: /content/work/cc/bin) on MY OWN
# copy of fs-build's tree, through my own copy of c-base's jail2.sh stand-in (uid 65534). Positive, then negative.
# No root Lake anywhere: the copy's .lake stays 65534-owned and Lake only runs inside the jail as 65534.
set -uo pipefail
V=/content/work/v-fs-statement; SRC=/content/work/fs-build/B; DST=$V/B; P=$DST/lean/complex588; L=$V/cmp; mkdir -p $L $V/cc/fresh-bin
TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin
echo "VCMP START $(date -u)"
until grep -q "VFS END" $V/run.out 2>/dev/null; do sleep 10; done
[ -e $DST ] && { echo "$DST exists; stop"; exit 1; }
s=$(date +%s); cp -a $SRC $DST || { echo COPY FAILED; exit 2; }; echo "copy done in $(( $(date +%s)-s ))s; owner of .lake: $(stat -c %U $P/.lake) mathlib: $(stat -c %U $P/.lake/packages/mathlib)"
cp /content/work/fs-build/cc/nounix $V/cc/nounix; cp /content/work/fs-build/cc/fresh-bin/lake $V/cc/fresh-bin/lake
sed "s#/content/work/fs-build/cc/nounix#$V/cc/nounix#" /content/work/fs-build/cc/jail2.sh > $V/cc/jail2.sh; chmod 755 $V/cc/jail2.sh $V/cc/nounix $V/cc/fresh-bin/lake
diff /content/work/fs-build/cc/jail2.sh $V/cc/jail2.sh
echo "tools: $(sha256sum /content/work/cc/bin/comparator /content/work/cc/bin/lean4export /content/work/cc/bin/nanoda_bin $V/cc/nounix $V/cc/fresh-bin/lake | awk '{print substr($1,1,16), $2}' | tr '\n' ' ')"
cd $P
echo "challenge sha256: $(sha256sum Work/FourierU8XMHH/UniformFourierChallenge.lean | cut -c1-64) neg: $(sha256sum Work/FourierU8XMHHNeg/UniformFourierChallenge.lean | cut -c1-64)"
for lab in vpos vneg; do
  C=comparator/fourier-u8xmhh.json; [ $lab = vneg ] && C=comparator/fourier-u8xmhh-neg.json
  mkdir -p $V/cc/scratch-$lab; chown -R 65534:65534 $V/cc/scratch-$lab
  w=0; while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 2 ] || [ "$(free -m | awk '/Mem:/ {print $7}')" -lt 14000 ]; do sleep 10; w=$((w+10)); done
  echo "=== $lab START $(date -u) (waited ${w}s) config $(sha256sum $C | cut -c1-16) $C"
  s=$(date +%s)
  env -i PATH=$TC:/content/work/cc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=$V/cc/home LANG=C.UTF-8 JAIL_TMP=$V/cc/scratch-$lab \
    /usr/bin/time -v -o $L/$lab.time $V/cc/jail2.sh $P /bin/sh -c "exec lake env /bin/sh -c 'PATH=$V/cc/fresh-bin:\$PATH; export PATH; exec comparator $C'" > $L/$lab.log 2>&1
  rc=$?
  echo "=== $lab END rc=$rc wall=$(( $(date +%s)-s ))s $(date -u)"
  grep -E "Building|accepts|okay|does not match|Const|rror|uses .sorry" $L/$lab.log | cut -c1-250 | tail -10
  grep -E "Maximum resident|Exit status" $L/$lab.time
done
echo "VCMP END $(date -u)"
