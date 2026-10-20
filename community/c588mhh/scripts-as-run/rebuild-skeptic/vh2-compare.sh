#!/bin/bash
# v-host-2: Shea's tools/Compare.lean (sha e48d6d3b...) on MY oleans: positive and negative; plus #print axioms.
B=/content/work/vh2/build
P=/content/work/mb/shea/lean/complex588/.lake/packages
MB=/content/work/mb/shea/lean/complex588/.lake/build/lib/lean
LEAN=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lean
L=/content/work/logs/v-host-2
C=/content/work/vh2/Compare.lean
cp /content/work/mb/shea/lean/complex588/tools/Compare.lean $C
sha256sum $C > $L/compare.sha
LP=$B:$MB
for p in mathlib batteries Qq aesop proofwidgets importGraph LeanSearchClient plausible Cli; do LP=$LP:$P/$p/.lake/build/lib/lean; done
export PATH=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin:$PATH
export LEAN_NUM_THREADS=2
cd /content/work/vh2
waitslot() { while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 2 ]; do sleep 5; done; }
waitslot
/usr/bin/time -v -o $L/cmp-pos.time env LEAN_PATH=$LP $LEAN --run $C Work.FourierC588MHH.UniformFourierChallenge Work.FourierC588MHH.Main OAI.PowerSaving.transform_mainC588MHH OAI.PowerSaving.convolution_mainC588MHH > $L/cmp-pos.log 2>&1; echo "pos rc=$?" >> $L/cmp.summary
waitslot
/usr/bin/time -v -o $L/cmp-neg.time env LEAN_PATH=$LP $LEAN --run $C Work.FourierC588MHHNeg.UniformFourierChallenge Work.FourierC588MHH.Main OAI.PowerSaving.transform_mainC588MHH OAI.PowerSaving.convolution_mainC588MHH > $L/cmp-neg.log 2>&1; echo "neg rc=$?" >> $L/cmp.summary
cat > /content/work/vh2/Ax.lean <<'EOA'
import Work.FourierC588MHH.Main
#check @OAI.PowerSaving.transform_mainC588MHH
#print OAI.PowerSaving.decimalExponentC588MHH
#print axioms OAI.PowerSaving.transform_mainC588MHH
#print axioms OAI.PowerSaving.convolution_mainC588MHH
#print axioms OAI.PowerSaving.XHost.engine_of_ginv_host
EOA
waitslot
/usr/bin/time -v -o $L/ax.time env LEAN_PATH=$LP $LEAN Ax.lean > $L/ax.log 2>&1; echo "ax rc=$?" >> $L/cmp.summary
echo DONE >> $L/cmp.summary
