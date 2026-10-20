#!/bin/bash
# h-L3-build: copy the .olean files my builds read into /dev/shm (RAM; ~2.1 GB) -- the Colab disk is I/O-saturated.
# Tree oleans come from the m-build tree (read only); my own Work/XHost and Work/FourierC588MH* stay on disk (symlinks).
set -e
S=/dev/shm/hl3; MB=/content/work/mb/shea/lean/complex588/.lake; MY=/content/work/hl3/shea/lean/complex588/.lake/build/lib/lean
mkdir -p $S/pk $S/tree
for p in mathlib batteries Qq aesop proofwidgets importGraph LeanSearchClient plausible Cli; do
  [ -d $MB/packages/$p/.lake/build/lib/lean ] || continue; mkdir -p $S/pk/$p
  (cd $MB/packages/$p/.lake/build/lib/lean && find . \( -name "*.olean*" -o -name "*.ir" -o -name "*.sig" \) -not -name "*.hash" -print0 | tar --null -cf - -T - | tar -xf - -C $S/pk/$p)
done
(cd $MB/build/lib/lean && find Work OAI WHTCheck \( -name "*.olean*" -o -name "*.ir" -o -name "*.sig" \) -not -name "*.hash" -print0 | tar --null -cf - -T - | tar -xf - -C $S/tree)
for d in XHost FourierC588MH FourierC588MHH FourierC588MHNeg FourierC588MHHNeg; do
  mkdir -p $MY/Work/$d; rm -rf $S/tree/Work/$d; ln -s $MY/Work/$d $S/tree/Work/$d
done
du -sh $S; touch $S/READY
