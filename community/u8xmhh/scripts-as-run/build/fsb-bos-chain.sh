#!/bin/bash
# fs-build: offline re-creation of the bigospeedrun/verifier check (verifier c92b9b6, profile "direct", challenge
# UniformFourier233.lean + Race.bound at the exact bound) for GXD.U8XMHH, with k-build's scripts and my paths.
# Starts only after one of k-build's two verifier chains has ended (so that at most 2 heavy Lean processes run), and
# each phase waits for fewer than 2 Lean processes (kb-bos-run's rule).   usage: fsb-bos-chain.sh
L=/content/work/logs/fs-build; F=/content/work/fs-build
until grep -q "^EXIT" /content/work/logs/k-build/bos-A.out 2>/dev/null || grep -q "^EXIT" /content/work/logs/k-build/bos-B.out 2>/dev/null; do sleep 30; done
echo "a k-build verifier chain has ended $(date -u): A: $(grep '^EXIT' /content/work/logs/k-build/bos-A.out) B: $(grep '^EXIT' /content/work/logs/k-build/bos-B.out)"
grep -q "Your solution is okay" $L/cmp/fsb-U8XMHH-pos.log || { echo "my positive comparator did not pass; verifier run skipped"; exit 1; }
bash $F/scripts/fsb-bos-prep.sh H $F/B/lean/complex588 Work.FourierU8XMHH.Main U8XMHH 16197705/10000000000 || { echo "PREP FAILED"; exit 2; }
bash $F/scripts/fsb-bos-run.sh H lean audit nanoda
