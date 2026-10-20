#!/bin/bash
# apply-overlay.sh <package dir> <tree = .../lean/complex588 of Ryan Shea's repository> [variant]
#
# Copies a package of package/ (package/u8x0hq, package/lever) into a clone of
# github.com/shea256/fourier-transform-below-nlogn at the commit named in the package's result.json.
# (Adapted from the script of the same name that I wrote for my earlier packages: the packages live in
# package/ instead of community/, and the lakefile lines that result.json lists under "lakefile_add" are appended
# to the tree's lakefile.lean when they are not there yet.)
#   * It first checks the package against its SHA256SUMS.
#   * mode "add" (result.json): every file must be NEW in the tree, or already equal to a file of this package
#     (another variant of it). If a file of the tree would be overwritten, nothing is copied.
#   * mode "replace": the package replaces files of the tree (c588-repack). Use a separate copy of the tree.
#   * variant: a folder under <package>/variants/ copied over the overlay after it (c588m: cyclic or packed).
#   * "requires" (result.json): packages that must be laid first (package/lever requires package/u8x0hq); each is
#     applied by this script, with its default files, before the package itself.
# Writes the list of copied files to $LOGDIR/applied.txt (default: ./logs/<package name>).
set -u
PKG=$(cd "${1:?package dir}" && pwd); T=$(cd "${2:?tree}" && pwd); V=${3:-}
NAME=$(basename "$PKG"); LOGDIR=${LOGDIR:-$PWD/logs/$NAME}; mkdir -p "$LOGDIR"
sha() { if command -v sha256sum >/dev/null; then sha256sum "$@"; else shasum -a 256 "$@"; fi; }
( cd "$PKG" && sha -c --quiet SHA256SUMS ) || { echo "SHA256SUMS mismatch: package damaged"; exit 91; }
[ -f "$T/lakefile.lean" ] && [ -d "$T/Work/Gauge" ] || { echo "not Shea's lean/complex588 tree: $T"; exit 92; }
for REQ in $(python3 -I -c 'import json,sys; print(" ".join(json.load(open(sys.argv[1])).get("requires",[])))' "$PKG/result.json"); do
  echo "== $NAME requires $REQ: applying package/$REQ first"
  LOGDIR="$LOGDIR/requires-$REQ" bash "$0" "$(dirname "$PKG")/$REQ" "$T" || { echo "could not apply $REQ"; exit 96; }
done
MODE=$(python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1])).get("mode","add"))' "$PKG/result.json")
WANT=$(python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["base"]["commit"])' "$PKG/result.json")
HEAD=$(git -C "$T" rev-parse HEAD 2>/dev/null || echo unknown)
[ "$HEAD" = "$WANT" ] || echo "WARNING: tree HEAD $HEAD, package built on $WANT"
if [ "$MODE" = add ]; then
  bad=0
  while read -r f; do
    d="$T/$f"
    if [ -e "$d" ] && ! cmp -s "$PKG/overlay/$f" "$d"; then
      ours=0
      for w in "$PKG"/variants/*/; do [ -e "$w$f" ] && cmp -s "$w$f" "$d" && ours=1; done
      [ $ours = 1 ] && continue
      echo "COLLISION: $f exists in the tree with other bytes"; bad=1
    fi
  done < <(cd "$PKG/overlay" && find . -type f | sed 's|^\./||' | sort)
  [ $bad = 0 ] || exit 93
fi
(cd "$PKG/overlay" && find . -type f | sed 's|^\./||' | sort) > "$LOGDIR/applied.txt"
cp -R "$PKG/overlay/." "$T/"
if [ -n "$V" ]; then
  [ -d "$PKG/variants/$V" ] || { echo "unknown variant $V"; exit 95; }
  cp -R "$PKG/variants/$V/." "$T/"
fi
python3 -I - "$PKG/result.json" "$T/lakefile.lean" "$LOGDIR/applied.txt" <<'PY'
import json, sys
add = json.load(open(sys.argv[1])).get("lakefile_add", [])
lf = open(sys.argv[2]).read(); have = [l.strip() for l in lf.splitlines()]
new = [l for l in add if l.strip() not in have]
if new:
    open(sys.argv[2], "a").write("\n" + "\n".join(new) + "\n")
with open(sys.argv[3], "a") as f:
    f.write("lakefile.lean: %d line(s) added: %s\n" % (len(new), new))
PY
echo "mode $MODE, variant ${V:-default}" >> "$LOGDIR/applied.txt"
echo "applied $(grep -c . "$LOGDIR/applied.txt") lines to $T; git status lines: $(git -C "$T" status --short --untracked-files=all | wc -l)"
