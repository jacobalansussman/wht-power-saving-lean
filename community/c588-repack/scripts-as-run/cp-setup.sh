#!/bin/bash
# c-pack: make the overlay tree /content/work/cpack (copy of the base tree) and put x-pack's B-quartic Lean data in it.
set -eu
B=/content/work/shea; T=/content/work/cpack; IN=/content/work/cpack-in
echo "SETUP START $(date -u)"
[ -d $T ] || { s=$(date +%s); cp -a $B $T; echo "copy took $(( $(date +%s)-s ))s"; }
mkdir -p $IN/x && tar -xzf $IN/cpack-in.tar.gz -C $IN/x
V=${1:-B-quartic}
( cd $IN/x/lean-out/$V && sha256sum -c logs/SHA256SUMS --quiet ) && echo "x-pack $V files: checksums OK ($(wc -l < $IN/x/lean-out/$V/logs/SHA256SUMS) files)"
P=$T/lean/complex588
cd $P
echo "base HEAD $(git rev-parse HEAD), status lines before overlay: $(git status --short | wc -l)"
cp -r $IN/x/lean-out/$V/Work/. Work/
git status --short
echo "changed/new files: $(git status --short | wc -l)"
cmp $IN/mine-12383291.lean Work/FourierC588/UniformFourierChallenge.lean && echo "challenge == my derivation from OpenAI fd4aeeb at 12383291"
cmp $IN/mine-12383292.lean Work/FourierC588Neg/UniformFourierChallenge.lean && echo "neg challenge == my derivation at 12383292"
grep -h "decimalExponentC588 : ℝ" Work/FourierC588/Goal.lean Work/FourierC588/UniformFourierChallenge.lean Work/FourierC588Neg/UniformFourierChallenge.lean
echo "SETUP END $(date -u)"
