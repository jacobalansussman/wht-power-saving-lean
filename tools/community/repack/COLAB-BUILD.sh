#!/usr/bin/env bash
# x-pack: Colab/Linux build of the tier-1 maximum on Ryan Shea's complex588 proof (Kothe's PR588 atm3 circuit).
# Upload x-pack-colab.tar.gz next to this script, then:   bash COLAB-BUILD.sh B-quartic     (or A-cubic)
#   B-quartic: optimal 5,009-bank table + new engine lemma Work.BlockAccounting.Moment4.term_le4
#              -> DFT saving 12383291/10^10 (kernel 12383292/10^10)
#   A-cubic:   optimal 5,009-bank table, unchanged engine (data swap only)
#              -> DFT saving 12383279/10^10 (kernel 12383280/10^10)
# Expected cost (Shea's own receipt for the same-size tables): fresh build ~35 min, peak RSS ~9 GB with LEAN_NUM_THREADS=1.
# This mirrors Shea's lean/complex588/scripts/build_and_check.sh WITHOUT its first two steps (check_provenance.py,
# regenerate.py): those assert byte-identity with HIS data and are expected to flag the replaced files by design.
set -euo pipefail
V="${1:-B-quartic}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="${WORK:-$HERE/xpack-build}"
mkdir -p "$WORK" && cd "$WORK"
[ -d shea ] || git clone https://github.com/shea256/fourier-transform-below-nlogn.git shea
# The Colab build used 5239484; its lean/ folder is identical to that of 0e61af6.
cd shea && git checkout -q 52394849c7ecf3412d39b925d75fc7fdb98994d8 && cd lean/complex588
command -v lake >/dev/null || { curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none; }
export PATH="$HOME/.elan/bin:$PATH"
mkdir -p /tmp/xpack && tar -xzf "$HERE/x-pack-colab.tar.gz" -C /tmp/xpack
(cd /tmp/xpack/lean-out/$V && sha256sum -c logs/SHA256SUMS --quiet) && echo "x-pack $V files: checksums OK"
cp -r /tmp/xpack/lean-out/$V/Work/. Work/
git status --short | head -60
LOG="$WORK/logs-$V"; mkdir -p "$LOG"
export LEAN_NUM_THREADS="${LEAN_NUM_THREADS:-1}"
lake exe cache get > "$LOG/cache.log" 2>&1
lake build Work.FourierC588.Main Work.FourierC588.UniformFourierChallenge \
  Work.FourierC588Neg.UniformFourierChallenge Work.GCert.Data.C588Controls \
  Work.FourierC588.Axioms > "$LOG/build.log" 2>&1 || { tail -40 "$LOG/build.log"; exit 1; }
lake env lean Work/FourierC588/Axioms.lean > "$LOG/axioms.log" 2>&1
lake env lean --run tools/Compare.lean Work.FourierC588.UniformFourierChallenge \
  Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 \
  OAI.PowerSaving.convolution_mainC588 > "$LOG/compare-fourier.log" 2>&1
if lake env lean --run tools/Compare.lean Work.FourierC588Neg.UniformFourierChallenge \
    Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 \
    OAI.PowerSaving.convolution_mainC588 > "$LOG/compare-negative.log" 2>&1; then
  echo "ERROR: the one-unit-stronger challenge was accepted"; exit 1
fi
python3 - "$LOG" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1])
assert 'RESULT: PASS' in (p/'compare-fourier.log').read_text()
bad=(p/'compare-negative.log').read_text()
assert 'definition body differs' in bad and 'decimalExponentC588' in bad
ax=(p/'axioms.log').read_text()
assert 'sorryAx' not in ax and ax.count('[propext, Classical.choice, Quot.sound]')==3
PY
grep -h "decimalExponentC588 : ℝ" Work/FourierC588/Goal.lean
echo "ALL CHECKS PASSED ($V)"
# Optional, the official comparator (as in Shea's RUNLOG; needs comparator + lean4export + nanoda built, landrun or shim):
#   lake env <comparator> comparator/fourier.json
