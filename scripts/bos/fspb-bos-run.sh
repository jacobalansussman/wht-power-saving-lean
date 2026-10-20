#!/bin/bash
# fsp-build: fs-build's fsb-bos-run.sh with my paths.  fs-build: k-build's kb-bos-run.sh with my paths (jail copy, stage, logs). k-build: m-build's mb-bos-run.sh (c-base's cb-bos-run.sh: the three phases of bigospeedrun/verifier verify():
# lean, audit, nanoda) on stage /content/work/kb/bos-<X>, through MY copy of jail2.sh.  The comparator compiles the
# challenge and every submission module itself (fresh .lake/build, one Lean job at a time).  Waits for a Lean slot
# (fewer than 2 lean processes on the machine) before each phase.   usage: kb-bos-run.sh <X> [phases]
X=$1; shift
B=/content/work/fsp-build/bos-$X; S=$B/stage; LOG=/content/work/logs/fsp-build/bos-$X; mkdir -p $LOG
TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin
PHASES=${*:-lean audit nanoda}
echo "threads for the solution build (lake adapter): $(cat $B/control/threads 2>/dev/null || echo 1)"
for ph in $PHASES; do
  ok=0; while [ $ok -lt 3 ]; do if [ "$(pgrep -c -x lean)" -ge 2 ]; then ok=0; sleep 15; else ok=$((ok+1)); sleep 4; fi; done
  if [ $ph != audit ]; then find $S/.lake/build/lib/lean $S/.lake/build/ir -maxdepth 1 -name 'VerifierChallenge.*' -print -delete 2>/dev/null; fi
  case $ph in
    lean|nanoda) CMD="exec lake env /bin/sh -c 'PATH=$B/control/bin:\$PATH; export PATH; exec comparator $B/control/$ph.json'";;
    audit) CMD="exec lake env lean VerifierAudit.lean";;
  esac
  s=$(date +%s)
  echo "=== $ph START $(date -u) epoch=$s olean files: $(find $S/.lake/build -name '*.olean' 2>/dev/null | wc -l)"
  echo "cmd: $CMD"
  env -i PATH=$TC:/content/work/cc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=/content/work/cc/home LANG=C.UTF-8 LEAN_NUM_THREADS=1 JAIL_TMP=$B/scratch \
    /usr/bin/time -v -o $LOG/$ph.time /content/work/fsp-build/cc/jail2.sh $S /bin/sh -c "$CMD" > $LOG/$ph.log 2>&1
  rc=$?; e=$(date +%s)
  echo "=== $ph END rc=$rc wall=$((e-s))s $(date -u) maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $LOG/$ph.time)"
  tail -4 $LOG/$ph.log
  echo "olean files after: $(find $S/.lake/build -name '*.olean' 2>/dev/null | wc -l)"
done
echo "BOS RUN END $(date -u)"
