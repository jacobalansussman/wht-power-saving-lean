#!/bin/bash
# fs-build (wht43 bound-check, 13 Oct 2026, Colab #2): my OWN tree for GXD.U8XMHH = k-combine's U8XM (rebuilt and
# checked on this Colab by k-build in /content/work/kb/B, Lake route, traces written) + h-L3-build's hosted network
# (Work/XHost, sources byte-identical to what h-L3-build compiled on Colab #1) + fs-price's generated U8XMHH data,
# glue and Fourier files.
#   1. waits until k-build's pipeline B has ended (its tree is then quiet), copies /content/work/kb/B (cp -a);
#   2. chowns MY copy's .lake back to root (kb's comparator left it to 65534; root Lake in a 65534-owned tree re-clones
#      and deletes packages), checks Lake sees the U8XM build as up to date and the packages intact;
#   3. checks the U8XM overlay in my copy against k-combine's package SHA256SUMS, unpacks the XHost and U8XMHH sources
#      and checks their hashes.
# usage: fsb-setup.sh     log: /content/work/logs/fs-build/setup.out
set -uo pipefail
. /content/work/env.sh
F=/content/work/fs-build; L=/content/work/logs/fs-build; mkdir -p $L
SRC=/content/work/kb/B; DST=$F/B; T=$DST/lean/complex588
PKG=/content/work/kb/pkg/B/k-combine-u8xm
echo "SETUP START $(date -u)"
until grep -qE "^EXIT|KB-RUN END" /content/work/logs/k-build/B.out 2>/dev/null; do sleep 20; done
echo "k-build B pipeline ended: $(grep -E '^EXIT|KB-RUN END' /content/work/logs/k-build/B.out | tr '\n' ' ') $(date -u)"
grep -E "accepts|okay|rc=" /content/work/logs/k-build/B.out | tail -6
[ -e $DST ] && { echo "$DST exists; stop"; exit 1; }
s=$(date +%s); cp -a $SRC $DST || { echo "COPY FAILED"; exit 2; }; echo "copy done in $(( $(date +%s)-s ))s; $(du -sh $DST | cut -f1)"
chown -R root:root $T/.lake && echo "my copy's .lake chowned to root"
cd $T
echo "HEAD $(git -C $DST rev-parse HEAD); tracked changes: $(git -C $DST status --porcelain --untracked-files=no | tr '\n' ' ')"
echo "mathlib before: $(git -C .lake/packages/mathlib rev-parse HEAD 2>&1) ; manifest: $(grep -A6 '"name": "mathlib"' lake-manifest.json | grep rev)"
ls .lake/packages
LEAN_NUM_THREADS=1 lake build --no-build Work.FourierU8XM.Main Work.FourierU8XM.UniformFourierChallenge Work.FourierU8XMNeg.UniformFourierChallenge > $L/pre-nobuild.log 2>&1; echo "lake --no-build U8XM rc=$? : $(tail -1 $L/pre-nobuild.log)"
echo "mathlib after: $(git -C .lake/packages/mathlib rev-parse HEAD 2>&1) ; packages: $(ls .lake/packages | tr '\n' ' ')"
# U8XM overlay files of my copy = k-combine's package
miss=0; n=0
while read h f; do case "$f" in overlay/*) n=$((n+1)); g=${f#overlay/}; h2=$(sha256sum $g | cut -d' ' -f1); [ "$h" = "$h2" ] || { echo "DIFFERS: $g"; miss=1; };; esac; done < $PKG/SHA256SUMS
echo "U8XM overlay files checked against k-combine's SHA256SUMS: $n, $([ $miss = 0 ] && echo all identical || echo MISMATCH)"
echo "k-combine package tarball: $(sha256sum /content/work/kb/pkg/k-combine-u8xm.tar.gz | cut -c1-64) (Mac: b18723496d576892e0bf907cc1f5a9d05335ac40ced9bb34b4199a49e6f3b5f9)"
[ -e Work/XHost ] && { echo "Work/XHost already present; stop"; exit 3; }
echo "48d9aaf1cd0763c932524a0243ca91607f61af33957df5f2ef5debb64bf43a4b  $F/in/xhost-src.tgz" | sha256sum -c - || exit 4
echo "80e197279e4c748b01ba0103aee213e76279132228ec81269b072fa1f517ef1b  $F/in/u8xmhh-src.tgz" | sha256sum -c - || exit 4
tar xzf $F/in/xhost-src.tgz && tar xzf $F/in/u8xmhh-src.tgz || exit 5
(sha256sum -c --quiet SHA256SUMS && echo "U8XMHH-SOURCES-OK ($(wc -l < SHA256SUMS) files)") || exit 6
sha256sum -c --quiet $F/in/xhost.sha256 && echo "XHOST-SOURCES-OK (= h-L3-build's)" || exit 7
mv SHA256SUMS $F/in/SHA256SUMS.u8xmhh; mv MODULES-U8XMHH.txt $F/in/
git -C $DST status --porcelain | head -5; echo "untracked entries: $(git -C $DST status --porcelain | wc -l)"
touch $F/SETUP-OK; echo "SETUP END $(date -u)"
