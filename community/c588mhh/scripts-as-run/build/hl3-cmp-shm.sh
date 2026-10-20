#!/bin/bash
# hl3-cmp-shm.sh (h-L3-build: m-build mb-cmp.sh / c-base cb-cmp.sh, log folder changed, and LEAN_PATH = the RAM-disk
# copies of the same .olean files (/dev/shm/hl3, see hl3-shm.sh) instead of `lake env`: the Colab disk is I/O-saturated) <label> <project dir> <config.json> [plain]
# One run of the official comparator (`lake env comparator <config>`, Shea's README), inside jail2.sh
# (stand-in for landrun/Landlock and systemd-run, which this kernel/container lack), timed with /usr/bin/time -v.
L=$1; T=$2; C=$3; MODE=${4:-jail}
LOG=/content/work/logs/h-L3-build/cmp
mkdir -p $LOG
TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin
LP=/dev/shm/hl3/tree
for p in mathlib batteries Qq aesop proofwidgets importGraph LeanSearchClient plausible Cli; do [ -d /dev/shm/hl3/pk/$p ] && LP=$LP:/dev/shm/hl3/pk/$p; done
cd "$T" || exit 90
echo "=== $L START $(date -u) epoch=$(date +%s) tree=$T config=$C mode=$MODE"
echo "git HEAD: $(git rev-parse HEAD)   git status lines: $(git status --short | wc -l)"
git status --short | head
echo "config: $(sha256sum "$C")"; cat "$C"
echo "tools: $(sha256sum /content/work/cc/bin/comparator /content/work/cc/bin/lean4export /content/work/cc/bin/landrun /content/work/cc/bin/nanoda_bin | tr "\n" " ")"
echo "olean files before the run: $(find .lake/build -name "*.olean" 2>/dev/null | wc -l)"
free -m | head -2
s=$(date +%s)
mkdir -p /content/work/cc/scratch-$L; chown 65534:65534 /content/work/cc/scratch-$L
if [ "$MODE" = plain ]; then W=""; else W="/content/work/cc/jail2.sh $T"; fi
env -i PATH=$TC:/content/work/cc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=/content/work/cc/home LANG=C.UTF-8 JAIL_TMP=/content/work/cc/scratch-$L \
  /usr/bin/time -v -o $LOG/$L.time $W env LEAN_PATH=$LP comparator "$C"
rc=$?
e=$(date +%s)
echo "=== $L END rc=$rc wall=$((e-s))s $(date -u) epoch=$e"
grep -E "Maximum resident|Elapsed|Exit status|User time|System time" $LOG/$L.time
echo "olean files after the run: $(find .lake/build -name "*.olean" 2>/dev/null | wc -l)"
echo "git status lines after: $(git status --short | wc -l)"
free -m | head -2
exit $rc
