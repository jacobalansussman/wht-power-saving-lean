#!/bin/bash
# i-skel: copy a staged overlay (Work/..., comparator/...) into my tree.  A file that already exists is SKIPPED if it is
# byte-identical and REFUSED (nothing copied) if it differs.  Writes $R/in-<TAG>.sha256 and $R/in-<TAG>.files
# (the .lean files, the list the sorry / forbidden greps read).
# usage: isk-install.sh <TAG> <staged dir>
set -uo pipefail
. ${ISK_SCRIPTS:-/content/work/i-skel/scripts}/isk-lib.sh
TAG=$1; IN=$2
cd $T || fail "no $T"
[ -d "$IN" ] || fail "no $IN"
bad=0; new=0; same=0
while read -r f; do
  if [ -e "$T/$f" ]; then
    if cmp -s "$IN/$f" "$T/$f"; then same=$((same+1)); else echo "REFUSE: $f exists and differs"; bad=1; fi
  else new=$((new+1)); fi
done < <(cd $IN && find . -type f | sed 's#^\./##' | sort)
[ $bad = 0 ] || fail "would overwrite a different file"
(cd $IN && find . -type f | sed 's#^\./##' | sort) | while read -r f; do
  [ -e "$T/$f" ] || { mkdir -p "$(dirname "$T/$f")"; cp -a "$IN/$f" "$T/$f"; }
done
chown -R root:root $T/Work $T/comparator 2>/dev/null
(cd $IN && find . -type f | sed 's#^\./##' | sort | xargs sha256sum) > $R/in-$TAG.sha256
(cd $IN && find . -type f -name '*.lean' | sed 's#^\./##' | sort) > $R/in-$TAG.files
echo "INSTALL $TAG: $new new, $same identical (skipped), from $IN; hashes $R/in-$TAG.sha256"
