#!/bin/bash
# run-comparator.sh <label> <tree = .../lean/complex588> <config.json> [plain]
#
# One run of the official comparator (`lake env comparator <config>`), as in Ryan Shea's README, timed with
# /usr/bin/time -v. On a machine without Landlock and systemd (the Colab VM) it runs inside the stand-in
# community/base-logs/scripts/jail2.sh: new mount, network and IPC namespaces, every mount read-only except
# <tree>/.lake and /dev, user nobody (65534), no capabilities, no_new_privs, and AF_UNIX sockets refused by
# the seccomp helper nounix (nounix.c). It is NOT landrun. With "plain" it runs without the stand-in.
#
# Environment: CC = folder with bin/{comparator,lean4export,landrun,nanoda_bin,nounix} and jail2.sh
# (built by community/base-logs/scripts/cb-tools.sh at the verifier's pins); TC = the bin folder of the
# Lean toolchain v4.34.1; LOGDIR (default ./logs/comparator).
#
# Ownership: the stand-in runs as nobody, so `chown -R 65534 <tree>/.lake` first. After that, do not run
# `lake` as root in that tree: git refuses the nobody-owned package clones ("dubious ownership"), Lake takes
# this for a changed URL, deletes .lake/packages and clones them again without their build outputs.
set -u
L=${1:?label}; T=$(cd "${2:?tree}" && pwd); C=${3:?config}; MODE=${4:-jail}
CC=${CC:-/content/work/cc}; TC=${TC:-/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin}
LOGDIR=${LOGDIR:-$PWD/logs/comparator}; mkdir -p "$LOGDIR"
cd "$T" || exit 90
echo "=== $L START $(date -u) tree=$T config=$C mode=$MODE"
echo "git HEAD: $(git rev-parse HEAD)   git status lines: $(git status --short | wc -l)"
echo "config: $(sha256sum "$C")"; cat "$C"
echo "tools: $(sha256sum $CC/bin/comparator $CC/bin/lean4export $CC/bin/nanoda_bin | tr "\n" " ")"
mkdir -p "$CC/scratch-$L"; chown 65534:65534 "$CC/scratch-$L"
if [ "$MODE" = plain ]; then W=""; else W="$CC/jail2.sh $T"; fi
s=$(date +%s)
env -i PATH=$TC:$CC/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=$CC/home LANG=C.UTF-8 JAIL_TMP=$CC/scratch-$L \
  /usr/bin/time -v -o "$LOGDIR/$L.time" $W lake env comparator "$C"
rc=$?
echo "=== $L END rc=$rc wall=$(( $(date +%s) - s ))s $(date -u)"
grep -E "Maximum resident|Elapsed|Exit status" "$LOGDIR/$L.time"
exit $rc
