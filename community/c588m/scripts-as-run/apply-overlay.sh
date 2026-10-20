#!/bin/bash
# apply-overlay.sh <tree = .../lean/complex588> [annealed|cyclic|packed]
# m-gen (13 Oct 2026).  Verifies the package (SHA256SUMS), then copies overlay/ into Ryan Shea's tree
# (shea256/fourier-transform-below-nlogn @ 5239484, lean/complex588).  Only NEW files are added: if a destination
# exists with different bytes the script stops before copying anything.  `cyclic` (fallback) or `packed` (optional
# upgrade) then copies variants/<v>/ over the overlay files it replaces (all of them files of this package);
# switching between the three variants is allowed.  Records what it did in
# /content/work/logs/${LBL:-m-c588m}/applied.txt.
set -u
T=${1:?tree}; V=${2:-annealed}
PKG=$(cd "$(dirname "$0")/.." && pwd)
L=/content/work/logs/${LBL:-m-c588m}; mkdir -p "$L"
cd "$PKG" && sha256sum -c --quiet SHA256SUMS || { echo "SHA256SUMS mismatch: package damaged"; exit 91; }
[ -f "$T/lakefile.lean" ] && [ -d "$T/Work/Gauge" ] || { echo "not Shea's lean/complex588 tree: $T"; exit 92; }
bad=0
while read -r f; do
  d="$T/$f"
  if [ -e "$d" ] && ! cmp -s "$PKG/overlay/$f" "$d"; then
    ours=0; for w in cyclic packed; do [ -e "$PKG/variants/$w/$f" ] && cmp -s "$PKG/variants/$w/$f" "$d" && ours=1; done
    [ $ours = 1 ] && continue                                                                         # our own other variant
    echo "COLLISION: $f exists in the tree with other bytes"; bad=1
  fi
done < <(cd "$PKG/overlay" && find . -type f | sed 's|^\./||' | sort)
[ $bad = 0 ] || exit 93
(cd "$PKG/overlay" && find . -type f | sed 's|^\./||' | sort) > "$L/applied.txt"
cp -R "$PKG/overlay/." "$T/"
case "$V" in
  annealed) ;;
  cyclic|packed) cp -R "$PKG/variants/$V/." "$T/" ;;
  *) echo "unknown variant $V"; exit 95 ;;
esac
echo "variant $V" >> "$L/applied.txt"
( cd "$T" && echo "git status lines (new files are untracked): $(git status --short --untracked-files=all | wc -l)" ) | tee -a "$L/applied.txt"
echo "applied $(grep -c . "$L/applied.txt") lines to $T ($V)"
