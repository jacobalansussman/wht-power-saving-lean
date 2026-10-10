#!/usr/bin/env python3
"""Fetch the nine outside input files of the rebuild, read-only, and check the sha256 of each.

usage: python3 -I -B tools/rebuild/fetch_inputs.py <directory> [pr233]

The files are data files (JSON) of github.com/CrocSwap/integer-mult-bounds at two fixed commits; rbpaths.py
lists them with their sizes and sha256.  For each file that is not yet in <directory> this script sends one
plain HTTPS GET to
    https://raw.githubusercontent.com/CrocSwap/integer-mult-bounds/<commit>/<path>
and writes the answer to <directory>/pr168/<path> or <directory>/pr193/<path>, but only if its sha256 is the
expected one.  A file that is already there is not fetched again; it is checked like the others.  So the
files may also be put there by other means (for example from a git checkout of the two commits) and this
script then only checks them.

With the word pr233 after the directory, one more file is fetched and checked in the same way, to
<directory>/pr233/<path>: the optional input of rbpaths.py, another choice of the hand-over pairs.  The rebuild
of the published certificate does not need it; rebuild.py reads it only when it is named with pairs=<file>.

Nothing is sent except the nine GET requests (ten with pr233): no credentials, no git, no sign-in.  Nothing
that is fetched is executed or imported.  Exit code 0 = all nine files (and, with pr233, the optional one) are
in <directory> with the expected sha256.
Python 3, standard library only.
"""
import hashlib
import os
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rbpaths

RAW = "https://raw.githubusercontent.com/%s/%s/%s"


def main():
    if len(sys.argv) not in (2, 3) or sys.argv[2:] not in ([], ["pr233"]):
        sys.exit(__doc__)
    dest = rbpaths.slash(sys.argv[1])
    extra = rbpaths.OPTIONAL if sys.argv[2:] else []
    for group, path, size, want, _ in rbpaths.FILES + extra:
        target = dest + group + "/" + path
        if os.path.exists(target):
            continue
        url = RAW % (rbpaths.HUB, rbpaths.COMMIT[group], path)
        print("GET %s" % url, flush=True)
        data = b""
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                while len(data) <= size:             # never more than one byte beyond the expected size
                    block = r.read(size + 1 - len(data))
                    if not block:
                        break
                    data += block
        except Exception as e:                       # no file is written; the check below reports it as MISSING
            print("  not fetched: %s" % e)
            continue
        got = hashlib.sha256(data).hexdigest()
        if got != want:
            print("  NOT WRITTEN: %d bytes with sha256 %s, expected %d bytes with %s" % (len(data), got, size, want))
            continue
        os.makedirs(os.path.dirname(target), exist_ok=True)
        with open(target, "wb") as f:
            f.write(data)
    print("input files in %s" % dest)
    bad = rbpaths.verify(dest)
    print("RESULT: %s" % ("all %d input files have the expected sha256" % len(rbpaths.FILES) if bad == 0
                          else "%d of %d input files are missing or different" % (bad, len(rbpaths.FILES))))
    if extra:
        print("optional input file in %s" % dest)
        bad2 = rbpaths.verify(dest, files=extra)
        print("RESULT: %s" % ("the optional input file has the expected sha256" if bad2 == 0
                              else "the optional input file is missing or different"))
        bad += bad2
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
