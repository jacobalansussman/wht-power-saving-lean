#!/bin/bash
# k-lower: the OFFICIAL comparator (c-base's tools: comparator d03acab, lean4export, nanoda; jail2.sh stand-in sandbox) on
# the U8X result.  Copy of the k-lower tree (U8X modules compiled with `lean -o`, NO Lake traces), so the comparator's
# own `lake build` (inside the jail, uid 65534, one job at a time through c-base's lake adapter, as c-base's run A3)
# recompiles the 73 new modules of the closure; then the comparator checks.  Then the negative config (one unit
# stronger, must fail).  Waits until fewer than 2 heavy Lean processes run before each run.
set -x
SRC=/content/work/k-lower; T=/content/work/cc/t-klower; P=$T/lean/complex588
LOGD=/content/work/logs/k-lower-u8x/cmp; mkdir -p $LOGD
s=$(date +%s)
[ -d $T ] || cp -a $SRC $T || exit 91
echo "copy wall=$(( $(date +%s)-s ))s"
cd $P
cat > comparator/fourier-u8x-neg.json <<JSON
{
  "challenge_module": "Work.FourierU8XNeg.UniformFourierChallenge",
  "solution_module": "Work.FourierU8X.Main",
  "theorem_names": ["OAI.PowerSaving.transform_mainU8X", "OAI.PowerSaving.convolution_mainU8X"],
  "definition_names": [],
  "permitted_axioms": ["propext", "Quot.sound", "Classical.choice"],
  "enable_nanoda": true
}
JSON
TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin
for L in klower-U8X-pos klower-U8X-neg; do
  C=comparator/fourier-u8x.json; [ $L = klower-U8X-neg ] && C=comparator/fourier-u8x-neg.json
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
