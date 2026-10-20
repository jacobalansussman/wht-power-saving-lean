#!/bin/sh
# check.sh [PROJECT] [UNIFORMFOURIER]
#
# Elaborates the two Lean files of this folder and prints their axioms.
#   PROJECT         the Lean project to elaborate in (default: the repository root, three levels up from
#                   research/existence-routes/programs/lean-existence-counts/); it must have built
#                   Work.GFrame.Engine.Main and Work.GFrame.Engine.Cert (`lake build Work.GFrame.Engine.Main
#                   Work.GFrame.Engine.Cert`) and Mathlib.
#   UNIFORMFOURIER  OpenAI's ComparatorChallenges/UniformFourier.lean (default: the copy in the repository's
#                   third-party/openai-math/lean/). Its sha256 must be the one below; ExLegitClosure.lean is
#                   that file followed by ExLegitClosure.tail, unchanged.
# Writes only into a fresh temporary folder. Expected: no errors; the only `sorry` warnings are the two in
# OpenAI's own file (its open challenge theorems, lines 326 and 330 of ExLegitClosure.lean); every
# `#print axioms` line lists only propext, Classical.choice and Quot.sound.
set -eu
HERE=$(cd "$(dirname "$0")" && pwd -P)
PROJECT=${1:-$HERE/../../../..}
UF=${2:-$PROJECT/third-party/openai-math/lean/ComparatorChallenges/UniformFourier.lean}
WANT=4c571f275bd7506c3ccd1e6d78f2fca77e07112c969c1c62ad8359d61544ec74
if command -v sha256sum >/dev/null 2>&1; then GOT=$(sha256sum "$UF" | cut -d' ' -f1); else GOT=$(shasum -a 256 "$UF" | cut -d' ' -f1); fi
[ "$GOT" = "$WANT" ] || { echo "UniformFourier.lean has sha256 $GOT, expected $WANT" >&2; exit 1; }
TMP=$(mktemp -d)
cat "$UF" "$HERE/ExLegitClosure.tail" > "$TMP/ExLegitClosure.lean"
cp "$HERE/ExLegit.lean" "$TMP/ExLegit.lean"
cd "$PROJECT"
echo "== ExLegitClosure.lean (imports Mathlib only)"
lake env lean "$TMP/ExLegitClosure.lean"
echo "== ExLegit.lean (imports the engine)"
lake env lean "$TMP/ExLegit.lean"
rm -rf "$TMP"
echo "== done"
