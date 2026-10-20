#!/bin/bash
# m-build: c-base's cb-bos-run.sh (the three phases of bigospeedrun/verifier verify(): lean, audit, nanoda) on the
# C588M stage made by mb-bos-prep.sh. Each phase in jail2.sh (stand-in; the verifier uses bwrap + systemd-run,
# unavailable here), env -i, /usr/bin/time -v. The comparator compiles the challenge and every submission module
# itself (fresh .lake/build). Waits while 2 or more Lean processes run (Colab rule).
B=/content/work/mb/bos; S=$B/stage; LOG=/content/work/logs/m-build/bos; mkdir -p $LOG
TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin
PHASES=${*:-lean audit nanoda}
echo "threads for the solution build (lake adapter): $(cat $B/control/threads 2>/dev/null || echo 1)"
/content/work/cb-sampler.sh $LOG/mem-sampler.log & SPID=$!
for ph in $PHASES; do
  ok=0; [ "${NOWAIT:-0}" = 1 ] && ok=3   # NOWAIT=1: phases run back to back on the Lean slot handed over by mb-final.sh
  while [ $ok -lt 3 ]; do if [ "$(pgrep -c lean)" -ge 2 ]; then ok=0; sleep 15; else ok=$((ok+1)); sleep 4; fi; done
  if [ $ph != audit ]; then find $S/.lake/build/lib/lean $S/.lake/build/ir -maxdepth 1 -name 'VerifierChallenge.*' -print -delete 2>/dev/null; fi
  case $ph in
    lean|nanoda) CMD="exec lake env /bin/sh -c 'PATH=$B/control/bin:\$PATH; export PATH; exec comparator $B/control/$ph.json'";;
    audit) CMD="exec lake env lean VerifierAudit.lean";;
  esac
  s=$(date +%s)
  echo "=== $ph START $(date -u) epoch=$s olean files: $(find $S/.lake/build -name '*.olean' 2>/dev/null | wc -l)"
  echo "cmd: $CMD"
  env -i PATH=$TC:/content/work/cc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=/content/work/cc/home LANG=C.UTF-8 LEAN_NUM_THREADS=1 JAIL_TMP=$B/scratch \
    /usr/bin/time -v -o $LOG/$ph.time /content/work/cc/jail2.sh $S /bin/sh -c "$CMD" > $LOG/$ph.log 2>&1
  rc=$?; e=$(date +%s)
  echo "=== $ph END rc=$rc wall=$((e-s))s $(date -u) maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $LOG/$ph.time)"
  tail -4 $LOG/$ph.log
  echo "olean files after: $(find $S/.lake/build -name '*.olean' 2>/dev/null | wc -l)"
done
kill $SPID
