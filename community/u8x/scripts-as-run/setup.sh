#!/bin/bash
# k-lower: own copy of the built tree (c-xexit = Shea's built lean/complex588 + x-exit-lean's XExit modules, built by
# colab-build), then the overlay (new files only; refuses to overwrite a different existing file).
set -euo pipefail
PKG=/content/work/k-lower-pkg-u8x
DST=/content/work/k-lower
SRC=${SRC:-/content/work/c-xexit}
echo "SETUP START $(date -u) src=$SRC"
cd $PKG && sha256sum -c --quiet SHA256SUMS.overlay && echo PACKAGE-OK
if [ ! -d $DST ]; then cp -a $SRC $DST; echo "copied $(du -sh $DST | cut -f1)"; fi
T=$DST/lean/complex588
cd $PKG/overlay
find . -type f | while read f; do
  d="$T/${f#./}"
  if [ -e "$d" ] && ! cmp -s "$f" "$d"; then echo "REFUSE: $d exists with other bytes"; exit 3; fi
  mkdir -p "$(dirname "$d")"; cp "$f" "$d"
done
grep -q 'lean_lib XExit' $T/lakefile.lean && echo "lakefile has XExit" || { echo "NO XExit lib"; exit 4; }
echo "SETUP END $(date -u)"
