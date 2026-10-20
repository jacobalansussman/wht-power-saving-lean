#!/bin/bash
# run-comparator.sh <label> <tree = .../lean/complex588> <config.json, relative to the tree> [plain]
#
# One run of the official comparator (`comparator <config>` inside `lake env`), timed with /usr/bin/time -v.
# (Adapted from the script of the same name that I wrote for my earlier packages, so that it runs exactly the
# recipe used for package/u8x0hq on 14 October 2026, scripts/bos/fspb-cmp.sh: the stand-in jail, the
# single-thread lake wrapper fresh-bin/lake, and a private scratch folder outside the tools folder.)
# On a machine without Landlock and systemd (the Colab VM) it runs inside the stand-in scripts/cc-jail/jail2.sh:
# new mount, network and IPC namespaces, every mount read-only except <tree>/.lake and a private temporary
# folder, user nobody (65534), no capabilities, no_new_privs, and AF_UNIX sockets refused by the seccomp helper
# nounix (nounix.c). It is NOT landrun. With "plain" it runs without the stand-in.
#
# Environment:
#   CC      folder with bin/{comparator,lean4export,nanoda_bin} (built by scripts/base/cb-tools.sh at the
#           verifier's pins); default /content/work/cc. Only read.
#   JAIL    folder with jail2.sh, nounix and fresh-bin/lake (scripts/cc-jail/; jail2.sh names the path of nounix
#           on its own lines); default /content/work/cc-jail
#   TC      the bin folder of the Lean toolchain v4.34.1; default /content/work/elan/toolchains/...
#   LOGDIR  default ./logs/comparator; the scratch folder of the run is $LOGDIR/scratch-<label>
#   CHOME   HOME inside the run (default $LOGDIR/home)
#
# Ownership: the stand-in runs as nobody, so this script runs `chown -R 65534 <tree>/.lake`. After that, do not
# run `lake` as root in that tree: git refuses the nobody-owned package clones ("dubious ownership"), Lake takes
# this for a changed URL, deletes .lake/packages and clones them again without their build outputs.
set -u
L=${1:?label}; T=$(cd "${2:?tree}" && pwd); C=${3:?config}; MODE=${4:-jail}
CC=${CC:-/content/work/cc}; JAIL=${JAIL:-/content/work/cc-jail}
TC=${TC:-/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin}
LOGDIR=${LOGDIR:-$PWD/logs/comparator}; mkdir -p "$LOGDIR"; LOGDIR=$(cd "$LOGDIR" && pwd)
SCR=$LOGDIR/scratch-$L; CHOME=${CHOME:-$LOGDIR/home}; mkdir -p "$SCR" "$CHOME"
cd "$T" || exit 90
echo "=== $L START $(date -u) tree=$T config=$C mode=$MODE"
echo "git HEAD: $(git rev-parse HEAD 2>/dev/null)   git status lines: $(git status --short 2>/dev/null | wc -l)"
echo "config: $(sha256sum "$C")"; cat "$C"
echo "tools: $(sha256sum $CC/bin/comparator $CC/bin/lean4export $CC/bin/nanoda_bin $JAIL/jail2.sh $JAIL/nounix $JAIL/fresh-bin/lake 2>/dev/null | awk '{print substr($1,1,16), $2}' | tr '\n' ' ')"
echo "olean files before: $(find .lake/build -name '*.olean' | wc -l)"
chown -R 65534:65534 "$T/.lake" "$SCR" "$CHOME"
s=$(date +%s)
if [ "$MODE" = plain ]; then
  env -i PATH=$TC:$CC/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=$CHOME LANG=C.UTF-8 \
    /usr/bin/time -v -o "$LOGDIR/$L.time" lake env comparator "$C"
else
  env -i PATH=$TC:$CC/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=$CHOME LANG=C.UTF-8 JAIL_TMP=$SCR \
    /usr/bin/time -v -o "$LOGDIR/$L.time" "$JAIL/jail2.sh" "$T" /bin/sh -c \
    "exec lake env /bin/sh -c 'PATH=$JAIL/fresh-bin:\$PATH; export PATH; exec comparator $C'"
fi
rc=$?
echo "=== $L END rc=$rc wall=$(( $(date +%s) - s ))s $(date -u)"
grep -E "Maximum resident|Elapsed|Exit status" "$LOGDIR/$L.time"
echo "olean files after: $(find .lake/build -name '*.olean' | wc -l)"
exit $rc
