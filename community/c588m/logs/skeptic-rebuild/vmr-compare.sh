#!/bin/bash
# v-m-rebuild: Compare.lean (ours, sha e48d6d3b) positive + negative against MY fresh build (out2 farm), toolchain bin on PATH.
. /content/work/env.sh
TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin; export PATH=$TC:$PATH
W=/content/work/v-m-rebuild; L=/content/work/logs/v-m-rebuild; B=/content/work/shea/lean/complex588; O2=$W/out2
LP=$O2; for p in $B/.lake/packages/*/; do LP=$LP:${p}.lake/build/lib/lean; done
export LEAN_PATH=$LP
echo "which lean: $(which lean)  Compare sha: $(sha256sum /content/work/ours/Compare.lean | cut -c1-16)"
echo "Main.olean resolved from: $(readlink -f $O2/Work/FourierC588M/Main.olean) (regular file: $([ -L $O2/Work/FourierC588M/Main.olean ] && echo no || echo yes))"
for k in pos neg; do
  ch=Work.FourierC588M.UniformFourierChallenge; [ $k = neg ] && ch=Work.FourierC588MNeg.UniformFourierChallenge
  while [ "$(pgrep -c lean)" -ge 2 ]; do sleep 10; done
  s=$(date +%s); /usr/bin/time -v -o $L/C-$k.time lean --run /content/work/ours/Compare.lean $ch Work.FourierC588M.Main OAI.PowerSaving.transform_mainC588M OAI.PowerSaving.convolution_mainC588M > $L/C-$k.log 2>&1
  echo "compare-$k rc=$? wall=$(( $(date +%s)-s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/C-$k.time) result: $(grep RESULT $L/C-$k.log)"
done
echo DONE
