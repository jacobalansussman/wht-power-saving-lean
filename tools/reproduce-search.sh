#!/bin/bash
# Re-run every search and pricing program of this repository and compare the results with the recorded ones
# (results/heads and results/g0). Python 3.9 or newer, standard library only. No Lean. About 15 minutes and at most
# 1 GB of memory on a laptop (14 minutes and 0.97 GB on an M-series Mac).
#
# usage: bash tools/reproduce-search.sh OUTDIR GATES_LEAN XHOST_DIR
#   OUTDIR      a new, empty folder for the fresh outputs
#   GATES_LEAN  Work/Bridge401/Gates.lean of Ryan Shea's tree (lean/complex588 at commit 5239484), sha256 077bcf03...
#   XHOST_DIR   the folder Work/XHost of package/u8x0hq/tree (HStage.lean, HCert.lean, HostK.lean are read)
#
# The last line says "ALL EQUAL" when every fresh output equals the recorded one (paths inside the JSON files
# are compared by file name only, and timing files are skipped).
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd -P)          # the repository root
OUT=$1; GATES=$2; XH=$3
[ -e "$OUT" ] && { echo "refusing: $OUT exists"; exit 1; }
mkdir -p "$OUT/heads" "$OUT/g0"
T=$R/tools/lever
P=$R/third-party/CrocSwap-integer-mult-bounds/staged544-terminal-gap-fill/owner
W=$P/inputs
py() { python3 -I -B "$@"; }
step() { echo "== $*"; }

# 0. the community's files are the pinned ones (git blob hashes, as GitHub stores them)
step "pinned community files"
(cd "$R/third-party/CrocSwap-integer-mult-bounds/staged544-terminal-gap-fill" &&
  while read sha f; do [ "$(git hash-object "$f")" = "$sha" ] || { echo "MISMATCH $f"; exit 1; }; done < BLOBS.txt)
echo "   all $(wc -l < "$R/third-party/CrocSwap-integer-mult-bounds/staged544-terminal-gap-fill/BLOBS.txt") files match their pins"

# 1. the community's figures, reproduced to 24 digits (two independent formulas)
step "controls.py";  py "$T/controls.py" > "$OUT/heads/controls.txt"
# 2. the lever's prices: T74 controls of the earlier plan, then #3945, then the other new words
step "lever_price.py"; py "$T/lever_price.py" > "$OUT/heads/lever-price.txt"
# 3. the head test on the community's own fill, its negative control, and a recount with no chronology
step "heads_fill.py"; py "$T/heads_fill.py" "$P" > "$OUT/heads/heads-3945.json"
py "$T/heads_fill.py" "$P" --control-bridge-first > "$OUT/heads/heads-3945-control.json"
step "recount.py";   py "$T/recount.py" "$P/FILL.json.gz" 360 > "$OUT/heads/recount-3945.txt"
# 4. our lowering of the two words of #3945 (and of #3815, the fallback) from the community's native canonical words
#    (data/native: the output of #3945's own native/run_native.py, run unchanged), with its checks and 7 mutants each
mkdir -p "$OUT/gcert"
for w in base extra; do
  step "lower_b.py #3945 $w"
  py "$T/lower_b.py" "$R/data/native/$w/CANONICAL-WORD.json.gz" "$OUT/gcert/S3945-$w.gcert1b.json.gz" "$OUT/g0/lower-$w.json" \
     --word "$W/WORD-$w.json.gz" --with-copies "$OUT/gcert/S3945-$w-copies.gcert1b.json.gz" > "$OUT/g0/lower-$w.txt"
  py "$T/mutants.py" "$R/data/native/$w/CANONICAL-WORD.json.gz" "$W/WORD-$w.json.gz" "$OUT/g0/mutants-$w.json" > "$OUT/g0/mutants-$w.txt"
  step "lower_b.py #3815 $w"
  py "$T/lower_b.py" "$R/data/native/old-$w/CANONICAL-WORD.json.gz" "$OUT/gcert/S3815-$w.gcert1b.json.gz" "$OUT/g0/lower-3815-$w.json" \
     --word "$W/old-WORD-$w.json.gz" --with-copies "$OUT/gcert/S3815-$w-copies.gcert1b.json.gz" > "$OUT/g0/lower-3815-$w.txt"
  py "$T/mutants.py" "$R/data/native/old-$w/CANONICAL-WORD.json.gz" "$W/old-WORD-$w.json.gz" "$OUT/g0/mutants-3815-$w.json" > /dev/null
done
step "lowered words byte-identical to data/gcert"
(cd "$OUT/gcert" && for f in *.gz; do cmp -s "$f" "$R/data/gcert/$f" || { echo "DIFFERS $f"; exit 1; }; done; echo "   $(ls | wc -l) files identical")
# 5. the head test in OUR lowering and OUR network order (about 2.5 minutes and 0.85 GB each)
step "g0_staged.py #3945"
py "$T/g0_staged.py" "$R" "$P" "$P/FILL.json.gz" "$R/data/gcert/S3945-base.gcert1b.json.gz" "$R/data/gcert/S3945-extra.gcert1b.json.gz" \
   "$P/PROFILE.json" "$OUT/g0/g0-3945.json" --label "#3945" --gates "$GATES" --xhost "$XH" > /dev/null
step "g0_staged.py #3815"
py "$T/g0_staged.py" "$R" "$P" "$W/old-FILL.json.gz" "$R/data/gcert/S3815-base.gcert1b.json.gz" "$R/data/gcert/S3815-extra.gcert1b.json.gz" \
   "$W/old-PROFILE.json" "$OUT/g0/g0-3815.json" --label "#3815" --gates "$GATES" --xhost "$XH" > /dev/null
step "g0_staged.py on the T74 regression data (data/t74ctl)"
D=$R/data/t74ctl
py "$T/g0_staged.py" "$R" "$D" "$D/FILL.json.gz" "$D/gauge44.gcert1b.json.gz" "$D/extra.gcert1b.json.gz" "$D/PROFILE.json" \
   "$OUT/g0/g0-T74-control.json" --gates "$GATES" --xhost "$XH" --label "T74 control" > "$OUT/g0/g0-T74-control.txt"
# 6. b of our lowering's profile by the community's moment formula, and the lever with our head count
step "price_ours.py"
py "$T/price_ours.py" "$OUT/g0/g0-3945.json" 0.001730663236534549309146 0.001730663236534549309147 > "$OUT/g0/price-3945.json"
py "$T/price_ours.py" "$OUT/g0/g0-3815.json" 0.001730578437266246818740 0.001730578437266246818741 > "$OUT/g0/price-3815.json"
# 7. compare
step "compare with the recorded results"
py "$R/tools/compare_results.py" "$OUT" "$R/results"
