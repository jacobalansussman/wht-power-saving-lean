#!/bin/bash
# v-fs-rebuild (skeptic, 13/14 Oct 2026, Colab #1): my OWN fresh copy of k-combine's ORIGINAL U8XM tree on Colab #1
# (/content/work/k-combine/lean/complex588; NOT k-build's Colab #2 rebuild that fs-build used), then the XHost and
# U8XMHH sources from the package tarballs (hash-checked), built from scratch.  log: /content/work/logs/v-fs-rebuild/setup.out
set -uo pipefail
. /content/work/env.sh
D=/content/work/v-fs-rebuild; L=/content/work/logs/v-fs-rebuild; mkdir -p $L
SRC=/content/work/k-combine/lean/complex588; T=$D/T
echo "SETUP START $(date -u)"
[ -e $T ] && { echo "$T exists; stop"; exit 1; }
echo "lean processes now: $(ps -C lean -o pid= | wc -l)"
s=$(date +%s); cp -a $SRC $T || { echo COPY FAILED; exit 2; }; echo "copy done $(( $(date +%s)-s ))s $(du -sh $T | cut -f1)"
chown -R root:root $T/.lake && echo "my copy's .lake chowned to root (source tree untouched: $(stat -c %U $SRC/.lake))"
cd $T
[ -e Work/XHost ] && { echo "Work/XHost already present in copy; stop"; exit 3; }
[ -e Work/FourierU8XMHH ] && { echo "FourierU8XMHH present; stop"; exit 3; }
echo "mathlib: $(git -C .lake/packages/mathlib rev-parse HEAD 2>&1); packages: $(ls .lake/packages | tr '\n' ' ')"
# U8XM overlay = k-combine package
miss=0; n=0
while read h f; do case "$f" in overlay/*) n=$((n+1)); g=${f#overlay/}; h2=$(sha256sum $g 2>/dev/null | cut -d' ' -f1); [ "$h" = "$h2" ] || { echo "DIFFERS: $g"; miss=1; };; esac; done < $D/in/kcombine-u8xm.SHA256SUMS
echo "U8XM overlay files vs k-combine SHA256SUMS: $n checked, $([ $miss = 0 ] && echo all identical || echo MISMATCH)"
LEAN_NUM_THREADS=1 lake build --no-build Work.FourierU8XM.Main > $L/pre-nobuild.log 2>&1; echo "lake --no-build U8XM rc=$? : $(tail -1 $L/pre-nobuild.log)"
echo "mathlib after: $(git -C .lake/packages/mathlib rev-parse HEAD 2>&1); packages: $(ls .lake/packages | tr '\n' ' ')"
echo "48d9aaf1cd0763c932524a0243ca91607f61af33957df5f2ef5debb64bf43a4b  $D/in/xhost-src.tgz" | sha256sum -c - || exit 4
echo "80e197279e4c748b01ba0103aee213e76279132228ec81269b072fa1f517ef1b  $D/in/u8xmhh-src.tgz" | sha256sum -c - || exit 4
tar xzf $D/in/xhost-src.tgz && tar xzf $D/in/u8xmhh-src.tgz || exit 5
sha256sum -c --quiet SHA256SUMS && echo "U8XMHH sources OK ($(wc -l < SHA256SUMS))" || exit 6
mv SHA256SUMS $D/in/SHA256SUMS.u8xmhh; mv MODULES-U8XMHH.txt $D/in/
echo "XHost files: $(ls Work/XHost | tr '\n' ' ')"
echo "olean present for XHost? $(ls .lake/build/lib/lean/Work/XHost 2>/dev/null | wc -l)  FourierU8XMHH? $(ls .lake/build/lib/lean/Work/FourierU8XMHH 2>/dev/null | wc -l)"
touch $D/SETUP-OK; echo "SETUP END $(date -u)"
