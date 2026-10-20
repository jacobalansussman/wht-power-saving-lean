#!/bin/bash
# k-combine: the OFFICIAL comparator (c-base's tools: comparator d03acab, lean4export 076e8e5, nanoda 3a24072; jail2.sh
# stand-in for landrun/systemd-run, exactly as k-lower's cmp-u8x.sh and c-base's cb-cmp.sh) on certificate NAME in
# tree T (our own copy; its .lake is given to uid 65534 as c-base did).  The comparator's own `lake build` inside the
# jail recompiles modules without Lake traces (one job at a time through c-base's lake adapter).  Positive config,
# then the negative config (challenge one unit stronger, must be REJECTED).  Waits for < 2 heavy Lean processes.
# usage: cmp.sh <NAME> <tree lean/complex588> <log dir>
set -x
N=$1; P=$2; LOGD=$3; n=$(echo $N | tr A-Z a-z)
mkdir -p $LOGD
cd $P
TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin
for L in kcomb-$N-pos kcomb-$N-neg; do
  C=comparator/fourier-$n.json; [ $L = kcomb-$N-neg ] && C=comparator/fourier-$n-neg.json
  mkdir -p /content/work/cc/scratch-$L
  chown -R 65534:65534 $P/.lake /content/work/cc/scratch-$L
  while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 2 ]; do sleep 10; done
  echo "=== $L START $(date -u) config: $(sha256sum $C)"
  s=$(date +%s)
  env -i PATH=$TC:/content/work/cc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=/content/work/cc/home LANG=C.UTF-8 JAIL_TMP=/content/work/cc/scratch-$L \
    /usr/bin/time -v -o $LOGD/$L.time /content/work/cc/jail2.sh $P /bin/sh -c "exec lake env /bin/sh -c 'PATH=/content/work/cc/fresh-bin:\$PATH; export PATH; exec comparator $C'" > $LOGD/$L.log 2>&1
  rc=$?
  echo "=== $L END rc=$rc wall=$(( $(date +%s)-s ))s $(date -u)"
  tail -5 $LOGD/$L.log
  grep -E "Maximum resident|Elapsed|Exit status" $LOGD/$L.time
done
echo CMPDONE
