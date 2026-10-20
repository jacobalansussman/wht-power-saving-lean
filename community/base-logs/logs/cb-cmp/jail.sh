#!/bin/bash
# jail.sh (agent c-comparator) -- NOT landrun, NOT Landlock. A stand-in built from Linux namespaces,
# ordinary user ids and a seccomp filter, because this kernel has no Landlock and the machine has no systemd.
#   usage: jail.sh <project directory> <command> [args...]
# The command runs
#   * in new mount, network and IPC namespaces (the new network namespace has no interface but a down loopback);
#   * with every mount point read-only, except /dev (and what is under it) and <project>/.lake;
#     /proc is left as the container has it;
#   * as user 65534:65534 (nobody), no supplementary groups, empty capability bounding set, no_new_privs;
#   * with socket(AF_UNIX) forbidden (stand-in for systemd RestrictAddressFamilies=~AF_UNIX).
# The environment is passed through unchanged.
set -eu
if [ "${1:-}" != "--inner" ]; then
  P=$(readlink -f "$1"); shift
  [ -d "$P/.lake" ] || { echo "jail: $P/.lake missing" >&2; exit 97; }
  exec unshare --mount --net --ipc --propagation private -- /bin/bash "$(readlink -f "$0")" --inner "$P" "$@"
fi
shift; P=$1; shift
mount --bind "$P/.lake" "$P/.lake"
fail=0
while read -r mp; do
  case "$mp" in
    /dev|/dev/*|/proc|/proc/*|"$P/.lake") continue ;;
  esac
  mount -o remount,bind,ro "$mp" 2>/dev/null || { echo "jail: could not make read-only: $mp" >&2; fail=1; }
done < <(awk "{print \$5}" /proc/self/mountinfo)
[ "$fail" = 0 ] || { echo "jail: refusing to run" >&2; exit 98; }
cd "$P"
exec setpriv --reuid 65534 --regid 65534 --clear-groups --bounding-set -all --inh-caps -all --no-new-privs -- /content/work/cc/bin/nounix "$@"
