#!/bin/bash
# fsp-build: the OFFICIAL comparator (comparator, lean4export, nanoda in /content/work/cc/bin) on U8X0HQ in MY tree, positive
# then negative, through MY copy of the stand-in jail (fs-build's / fsp-cc's jail2.sh, only the nounix path changed:
# namespaces, read-only mounts except <project>/.lake and a private /tmp, uid 65534, no capabilities, no_new_privs,
# socket(AF_UNIX) blocked) -- NOT the official landrun/systemd sandbox.  fsp-cc's fspcc-cmp.sh with my paths.
# Runs only when Compare.lean gave PASS / FAIL.  After the chown, no root Lake in this tree.
# usage: fspb-cmp.sh U8X0HQ
F=/content/work/fsp-build; P=$F/T/lean/complex588; TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin
heavy() { ps -C lean -o rss= | awk '$1>2000000' | wc -l; }
# my jail copy
if [ ! -e $F/cc/jail2.sh ]; then
  mkdir -p $F/cc/fresh-bin && cp -a /content/work/fsp-cc/cc/nounix /content/work/fsp-cc/cc/nounix.c $F/cc/ && cp -a /content/work/fsp-cc/cc/fresh-bin/lake $F/cc/fresh-bin/
  sed 's#^\# fsp-cc copy (14 Oct 2026) of fs-build/k-build copy#\# fsp-build copy (14 Oct 2026) of fsp-cc/fs-build/k-build copy#; s#/content/work/fsp-cc/cc/nounix#/content/work/fsp-build/cc/nounix#' /content/work/fsp-cc/cc/jail2.sh > $F/cc/jail2.sh; chmod 755 $F/cc/jail2.sh
fi
echo "f90dbc3996f3c56ba0b1d0fcec6331fc25a32780c002f38ef130b7f40701374e  $F/cc/nounix.c" | sha256sum -c - || exit 2
echo "1e6528c3f2c8b10b04f5c5e5a77d4ef71f17a09d35cbe836b14cb3022f1c6779  $F/cc/nounix" | sha256sum -c - || exit 2
echo "01f56ebe57de64f1bea763882a8910783ec31336c1d4b8577d497d1a7c5777c7  $F/cc/fresh-bin/lake" | sha256sum -c - || exit 2
diff /content/work/fsp-cc/cc/jail2.sh $F/cc/jail2.sh
cd $P
echo "tools: $(sha256sum /content/work/cc/bin/comparator /content/work/cc/bin/lean4export /content/work/cc/bin/nanoda_bin $F/cc/jail2.sh $F/cc/nounix $F/cc/fresh-bin/lake | awk '{print substr($1,1,16), $2}' | tr '\n' ' ')"
for N in "$@"; do
  n=$(echo $N | tr A-Z a-z); LD=/content/work/logs/fsp-build/$n; LOGD=$LD/cmp; mkdir -p $LOGD
  grep -q "RESULT: PASS" $LD/compare-pos.log 2>/dev/null && grep -q "RESULT: FAIL" $LD/compare-neg.log 2>/dev/null || { echo "$N: Compare.lean results not PASS/FAIL: skipped"; continue; }
  for L in $n-pos $n-neg; do
    C=comparator/fourier-$n.json; [ $L = $n-neg ] && C=comparator/fourier-$n-neg.json
    mkdir -p $F/cc/scratch-$L
    chown -R 65534:65534 $P/.lake $F/cc/scratch-$L
    w=0; while [ "$(heavy)" -ge 2 ] || [ "$(free -m | awk '/Mem:/ {print $7}')" -lt 14000 ]; do sleep 10; w=$((w+10)); done
    echo "=== $L START $(date -u) (slot waited ${w}s) config: $(sha256sum $C)"; cat $C
    echo "olean files before: $(find .lake/build -name '*.olean' | wc -l)"
    s=$(date +%s)
    env -i PATH=$TC:/content/work/cc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=/content/work/cc/home LANG=C.UTF-8 JAIL_TMP=$F/cc/scratch-$L \
      /usr/bin/time -v -o $LOGD/$L.time $F/cc/jail2.sh $P /bin/sh -c "exec lake env /bin/sh -c 'PATH=$F/cc/fresh-bin:\$PATH; export PATH; exec comparator $C'" > $LOGD/$L.log 2>&1
    rc=$?
    echo "=== $L END rc=$rc wall=$(( $(date +%s)-s ))s $(date -u)"
    grep -E "Building|Built|accepts|okay|does not match|Const|rror" $LOGD/$L.log | cut -c1-300 | tail -10
    grep -E "Maximum resident|Elapsed|Exit status" $LOGD/$L.time
    echo "olean files after: $(find .lake/build -name '*.olean' | wc -l)"
  done
done
echo CMPDONE
