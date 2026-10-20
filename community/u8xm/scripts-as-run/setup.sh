#!/bin/bash
# k-combine: apply the overlay of NEW files to a built tree (refuses to overwrite a different existing file; identical
# files, e.g. k-lower's Work/KLow/*.lean already in a k-lower copy, are accepted).  The base tree must be Ryan Shea's
# BUILT lean/complex588 (5239484) with x-exit-lean's XExit library (lakefile line `lean_lib XExit ...`), e.g.
# colab-build's /content/work/c-xexit or a copy of /content/work/k-lower.
# usage: setup.sh <unpacked package dir> <tree: .../lean/complex588>
set -euo pipefail
PKG=$1; T=$2
echo "SETUP START $(date -u) pkg=$PKG tree=$T"
cd $PKG && sha256sum -c --quiet SHA256SUMS.overlay && echo PACKAGE-OK
grep -q 'lean_lib XExit' $T/lakefile.lean && echo "lakefile has XExit" || { echo "NO XExit lib in $T/lakefile.lean"; exit 4; }
cd $PKG/overlay
n=0; same=0
while read f; do
  d="$T/${f#./}"
  if [ -e "$d" ]; then
    if cmp -s "$f" "$d"; then same=$((same+1)); continue; fi
    echo "REFUSE: $d exists with other bytes"; exit 3
  fi
  mkdir -p "$(dirname "$d")"; cp "$f" "$d"; n=$((n+1))
done < <(find . -type f)
echo "copied $n new files, $same identical files already present"
echo "SETUP END $(date -u)"
