#!/bin/bash
# i-skel SCRIPT CONTROL: my parameterised scripts on the EXISTING, prebuilt U8X0HQ instance of ours/T (suffix U8X0HQ,
# DFT figure 1 - 16256629/10^10).  They must reproduce fsp-build's verdicts of 14 Oct: challenges byte-identical,
# Compare.lean positive PASS / negative FAIL, official comparator accepted / rejected, bigospeedrun re-creation
# (lean, audit, nanoda) pass; plus a negative bigospeedrun run at one grid step too high (rejected).
# usage: isk-control.sh [bos threads, default 3]
set -uo pipefail
. ${ISK_SCRIPTS:-/content/work/i-skel/scripts}/isk-lib.sh
THR=${1:-3}
echo "CONTROL START $(date -u) on $(hostname)"
bash $S/isk-tree.sh || fail tree
bash $S/isk-build.sh U8X0HQ Work.FourierU8X0HQ.Main Work.FourierU8X0HQ.UniformFourierChallenge Work.FourierU8X0HQNeg.UniformFourierChallenge Work.FourierU8X0HQ.Axioms || fail "U8X0HQ targets"
bash $S/isk-checks.sh U8X0HQ 16256629 || fail "checks U8X0HQ"
bash $S/isk-cmp.sh U8X0HQ || fail "comparator U8X0HQ"
bash $S/isk-bos.sh Q U8X0HQ 16256629 $THR pos || fail "bos Q"
bash $S/isk-bos.sh Qneg U8X0HQ 16256629 $THR neg || echo "bos negative did not give the expected verdict (see log)"
ls -la $R/*-OK
echo "CONTROL END $(date -u)"
echo EXIT 0
