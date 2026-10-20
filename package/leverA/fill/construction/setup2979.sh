#!/bin/sh
# o-model (wht45/reopt), on C1: lay out #2979's own allocation code beside the scratch tree S/ it expects.
# Hub code is copied unchanged from main 20ecf2891 (worktree /content/work/reopt-o-model/hub-main); only symlinks are added.
set -e
B=/content/work/reopt-o-model
H=$B/hub-main/research
C=$B/c2979
rm -rf $C; mkdir -p $C/code $C/S/r13/vendor/pr2917/research $C/S/r12/vendor/pr2857/research $C/S/r12/vendor/pr2833/research $C/runs
cp $H/rank7-joint-endframes/construction/allocation/*.py $C/code/
cp $H/rank7-joint-endframes/tools/relax.py $H/rank7-joint-endframes/tools/root.py $C/code/
ln -s $H/paired-helper-extension $C/S/r13/vendor/pr2917/research/paired-helper-extension
ln -s $H/shared-kernel-repacking $C/S/r12/vendor/pr2857/research/shared-kernel-repacking
ln -s $H/gmx-t74-borrow43-handoff $C/S/r12/vendor/pr2833/research/gmx-t74-borrow43-handoff
cd $C/code && sha256sum *.py
cd $H/rank7-joint-endframes && grep -E '"tools/(relax|root).py"' PINS.json
