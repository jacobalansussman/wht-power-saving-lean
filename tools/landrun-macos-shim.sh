#!/bin/sh
# Stand-in for `landrun` (Linux Landlock sandbox, unavailable on macOS) so the official
# leanprover/comparator can run here. It translates landrun's path flags into a macOS
# seatbelt profile: no network, writes only under the paths comparator marks writable
# (plus /dev and the temp dirs). Set LANDRUN_SHIM_NOSANDBOX=1 to run the command directly.
rw=""
while [ $# -gt 0 ]; do
  case "$1" in
    --) shift; break ;;
    --rw|--rwx) p=$(cd "$2" 2>/dev/null && pwd -P || echo "$2"); rw="$rw (subpath \"$p\")"; shift 2 ;;
    --ro|--rox|--env) shift 2 ;;
    *) shift ;;
  esac
done
if [ -n "$LANDRUN_SHIM_NOSANDBOX" ]; then exec "$@"; fi
profile="(version 1)(allow default)(deny network*)(deny file-write*)(allow file-write* (subpath \"/dev\") (subpath \"/private/tmp\") (subpath \"/private/var/folders\") $rw)"
exec /usr/bin/sandbox-exec -p "$profile" "$@"
