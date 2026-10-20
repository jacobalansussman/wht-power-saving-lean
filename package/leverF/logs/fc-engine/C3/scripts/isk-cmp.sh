#!/bin/bash
# i-skel: the OFFICIAL comparator (comparator + lean4export + nanoda, /content/work/cc/bin) on instance SUF in my tree,
# positive (comparator/fourier-<suf>.json) then negative (comparator/fourier-<suf>-neg.json), through prep's stand-in
# jail /content/work/cc-jail/jail2.sh (namespaces, read-only mounts except <project>/.lake and a private /tmp,
# uid 65534, no capabilities, no_new_privs, socket(AF_UNIX) blocked) -- NOT the official landrun/systemd sandbox.
# Recipe: fsp-build fspb-cmp.sh / e-r45m1 mxt-build.sh.  Runs only after isk-checks.sh wrote CHECKS-<SUF>-OK.
# My tree's .lake is chowned to 65534 for the run and back to root afterwards.
# Writes $R/CMP-<SUF>-OK when positive = "Your solution is okay!" and negative = rejected.
# usage: isk-cmp.sh <SUF>
set -uo pipefail
. ${ISK_SCRIPTS:-/content/work/i-skel/scripts}/isk-lib.sh
SUF=$1; n=$(lower $SUF)
L=$LOGROOT/$n/cmp; mkdir -p $L
rm -f $R/CMP-$SUF-OK
[ -e $R/CHECKS-$SUF-OK ] || fail "no CHECKS-$SUF-OK: Compare.lean verdicts first"
(cd $J && sha256sum -c --quiet HASHES) || fail "jail kit hashes"
[ -x /content/work/cc/bin/comparator ] || fail "no comparator"
echo "tools: $(sha256sum /content/work/cc/bin/comparator /content/work/cc/bin/lean4export /content/work/cc/bin/nanoda_bin $J/jail2.sh $J/nounix $J/fresh-bin/lake | awk '{print substr($1,1,16), $2}' | tr '\n' ' ')"
cd $T || fail "no $T"
ok=1
for LB in $n-pos $n-neg; do
  C=comparator/fourier-$n.json; [ $LB = $n-neg ] && C=comparator/fourier-$n-neg.json
  [ -e $C ] || fail "no $C"
  mkdir -p $R/cc/scratch-$LB $R/cc/home
  chown -R 65534:65534 $T/.lake $R/cc/scratch-$LB
  wl 14000
  echo "=== $LB START $(date -u) config: $(sha256sum $C | cut -c1-16)"; cat $C
  echo "olean files before: $(find .lake/build -name '*.olean' | wc -l)"
  s=$(date +%s)
  env -i PATH=$TC:/content/work/cc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=$R/cc/home LANG=C.UTF-8 JAIL_TMP=$R/cc/scratch-$LB \
    /usr/bin/time -v -o $L/$LB.time $J/jail2.sh $T /bin/sh -c "exec lake env /bin/sh -c 'PATH=$J/fresh-bin:\$PATH; export PATH; exec comparator $C'" > $L/$LB.log 2>&1
  rc=$?
  echo "=== $LB END rc=$rc wall=$(( $(date +%s)-s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$LB.time) $(date -u)"
  grep -E "accepts|okay|does not match|Const|rror" $L/$LB.log | cut -c1-300 | tail -8
  echo "olean files after: $(find .lake/build -name '*.olean' | wc -l)"
  if [ $LB = $n-pos ]; then
    [ $rc = 0 ] && grep -q "Your solution is okay!" $L/$LB.log && echo "VERDICT comparator positive: accepted" || { echo "VERDICT comparator positive: NOT accepted"; ok=0; }
  else
    [ $rc != 0 ] && ! grep -q "Your solution is okay!" $L/$LB.log && echo "VERDICT comparator negative: rejected (expected)" || { echo "VERDICT comparator negative: NOT rejected"; ok=0; }
  fi
done
chown -R root:root $T/.lake
[ $ok = 1 ] && { touch $R/CMP-$SUF-OK; echo "CMP $SUF ALL AS EXPECTED"; echo EXIT 0; } || { echo "CMP $SUF: SOMETHING UNEXPECTED"; echo EXIT 1; exit 1; }
