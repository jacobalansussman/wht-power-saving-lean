#!/bin/bash
# m-build (13 Oct 2026): build m-gen's C588M package (copy-free #588) ON TOP of c-base's built copy of
# Ryan Shea's tree, in a private copy (/content/work/mb/shea), without rebuilding Mathlib or Shea's modules.
# usage (inside a fresh tmux server): mb-build.sh <base repo dir, e.g. /content/work/shea> [annealed|cyclic|packed]
# Logs: /content/work/logs/m-build/ (one log + /usr/bin/time -v per step).
set -u
. /content/work/env.sh
BASE=${1:-/content/work/shea}; V=${2:-annealed}
W=/content/work/mb; LBL=${LBL:-m-build}; L=/content/work/logs/$LBL; mkdir -p $W $L
R=$W/shea; T=$R/lean/complex588
PKG=/content/work/m-gen/c588m
export LEAN_NUM_THREADS=${LEAN_NUM_THREADS:-2}
T0=$(date +%s)
echo "MB START $(date -u) epoch=$T0 base=$BASE variant=$V"
/content/work/cb-sampler.sh $L/mem-sampler.log & SPID=$!
trap 'kill $SPID 2>/dev/null' EXIT
if [ ! -d $R ]; then s=$(date +%s); cp -a $BASE $R || exit 91; echo "copy wall=$(( $(date +%s)-s ))s"; fi
echo "copy HEAD $(git -C $R rev-parse HEAD); tracked changes: $(git -C $R status --short --untracked-files=no | wc -l); olean files: $(find $T/.lake/build -name '*.olean' | wc -l)"
(cd $PKG && sha256sum -c --quiet SHA256SUMS && echo PACKAGE-OK) || exit 89
LBL=$LBL bash $PKG/scripts/apply-overlay.sh $T $V || exit 92
bash /content/work/mb/build-c588m-mb.sh $T $LBL; BRC=$?
echo "build-c588m rc=$BRC wall-so-far=$(( $(date +%s)-T0 ))s"
cd $T || exit 90
step() { # step <name> <command...>
  local n=$1; shift; local s=$(date +%s)
  echo "=== STEP $n START $(date -u) cmd: $*"
  /usr/bin/time -v -o $L/$n.time "$@" > $L/$n.log 2>&1; local rc=$?
  echo "=== STEP $n END rc=$rc wall=$(( $(date +%s) - s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$n.time) last: $(tail -1 $L/$n.log)"
  return $rc
}
waitlean() { local ok=0; while [ $ok -lt 3 ]; do if [ "$(pgrep -c lean)" -ge 2 ]; then ok=0; sleep 15; else ok=$((ok+1)); sleep 4; fi; done; }   # m-build: 3 consecutive samples below 2 (gaps between other agents' modules are not free slots)
if lake build --no-build Work.FourierC588M.Main >/dev/null 2>&1; then
  waitlean
  step ours-compare-fourier lake env lean --run /content/work/ours/Compare.lean Work.FourierC588M.UniformFourierChallenge \
    Work.FourierC588M.Main OAI.PowerSaving.transform_mainC588M OAI.PowerSaving.convolution_mainC588M
  waitlean
  step ours-compare-negative lake env lean --run /content/work/ours/Compare.lean Work.FourierC588MNeg.UniformFourierChallenge \
    Work.FourierC588M.Main OAI.PowerSaving.transform_mainC588M OAI.PowerSaving.convolution_mainC588M
  echo "(negative control must be nonzero)"
  step nobuild lake build --no-build Work.FourierC588M.Main Work.FourierC588M.UniformFourierChallenge \
    Work.FourierC588MNeg.UniformFourierChallenge Work.GCert.Data.C588MControls Work.FourierC588M.Axioms
fi
# source scan: every file of the overlay as it sits in the tree
{
  echo "# grep of the overlay files in the tree for sorry/admit/axiom/native_decide/implemented_by/extern/unsafe/opaque (word match, comments included)"
  (cd $PKG/overlay && find . -name '*.lean' | sed 's|^\./||' | sort) | while read -r f; do
    grep -nwE 'sorry|admit|axiom|native_decide|implemented_by|extern|unsafe|opaque|sorryAx' "$T/$f" | sed "s|^|$f:|"
  done
  echo "# declarations starting with 'axiom' anywhere under OAI Work WHTCheck:"
  grep -rnE '^\s*(private |protected )?axiom ' OAI Work WHTCheck | head -20
  echo "# files of the tree changed vs git HEAD (tracked):"
  git status --short --untracked-files=no
  echo "# untracked files: $(git status --short --untracked-files=all | grep -c '^??')"
} > $L/source-scan.txt 2>&1
echo "olean count: $(find .lake/build/lib/lean -name '*.olean' | wc -l)"
echo "MB END $(date -u) total_wall=$(( $(date +%s)-T0 ))s build_rc=$BRC"
awk '{ if ($5+0>m) m=$5+0; if ($4+0>s) s=$4+0 } END { print "sampler: max single process RSS " m " MB, max summed lean/lake RSS " s " MB (all agents)" }' $L/mem-sampler.log
exit $BRC
