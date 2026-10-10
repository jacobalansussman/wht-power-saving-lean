#!/usr/bin/env python3
"""Rebuild the certificate tools/certificate/gcert1-p11-pr193.json.gz from the outside data files.

PRINCIPLE.  No program of the outside repository (github.com/CrocSwap/integer-mult-bounds) is executed or
imported, by this script or by any script it starts.  The outside files are JSON.  Their sha256 is checked
first, and then they are loaded as data with json.load.  The outside scripts were read as text.  Six scripts
of this directory are adapted from them (compile.py, graphgen.py, aligned.py, stage.py, and parts of build.py
and f2.py); each says so in its first lines, and NOTICE, section 6, names the outside files and their authors.

usage (from anywhere):
    python3 -I -B tools/rebuild/fetch_inputs.py <inputs>         once: fetches and checks the nine outside files
                                                                 (with the word pr233 added: also the optional tenth)
    python3 -I -B tools/rebuild/rebuild.py <inputs> [out=<file.json.gz>] [--keep <dir>] [pairs=<file> name=<name>]

  <inputs>       the directory that fetch_inputs.py filled
  out=<file>     also copy the rebuilt certificate to <file>, which must not exist yet
  --keep <dir>   work in <dir> (empty or absent) and keep it; default: a temporary directory, removed at the end
  pairs=<file> name=<name>
                 use another file of hand-over pairs (same format: {"pairs": [[donor, recipient, operation
                 index or null], ...]} on the aligned word) in place of the file of #193.  <name> (lower-case
                 letters and digits, not pr193) replaces pr193 in the names of the steps and goes into the
                 certificate.  There is then no published value to compare with: exit code 0 means that
                 gx.check1 accepted the result in step 6.  The file is checked for form and range in step 2
                 (aligned.py, load_pairs).  Whether its hand-overs are legal is decided by the frame scan of
                 step 4 and by gx.check1 in step 6: an illegal file stops the rebuild there, and no
                 certificate is written.  The nine standard files are checked as always.  rbpaths.py lists
                 one such file, of pull request #233, which fetch_inputs.py <inputs> pr233 puts at
                 <inputs>/pr233/research/source-assisted-v4-layer/witness/physical-pairs.json.

Steps, each a separate process (nothing is written into the source tree):
    0  rbpaths.verify        sha256 of the nine input files; nothing is run if one is missing or different
    1  build.py pr168        the helper word of pull request #168 from the modules and the carrier arcs
    2  aligned.py pr193      the aligned word of #193: arcs and frames transported, hand-over pairs loaded
    3  rewrite.py pr193      the in-place pair steps: one explicit gate list, named pr193xword
    4  phys.py pr193xword    stage word, frame scan of every register, block counts, number of unit moves
    5  emit.py pr193xword    the circuit as a file of format gcert/0
    6  tools/gx/gxconv.py    gcert/0 -> gcert/1; the result is tested by gx.check1 before it is written
The sha256 of the uncompressed JSON is then compared with the value in tools/gx/ORIGIN.md.  The sha256 of the
.gz file itself differs from run to run, because gzip stores a time stamp.
Exit code 0 = the uncompressed content is that of the published certificate.

The name pr193xword is part of the certificate (gxconv.py copies it into the field `derived_from`), so the
names of the steps must stay as they are if the content is to come out the same.

Needs Python 3.9 at least (tools/gx/gx.py imports math.lcm); standard library only.  No Lean is run.
"""
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rbpaths

CERT = "gcert1-p11-pr193.json.gz"
PUBLISHED = "tools/certificate/" + CERT
EXPECTED = "c64d2a95431c3cf24452de7391cce780bd71521da44ce23923ebb235f6c7962c"      # sha256 of the uncompressed JSON


def sha_unpacked(path):
    h, n = hashlib.sha256(), 0
    with gzip.open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
            n += len(block)
    return h.hexdigest(), n


def main():
    args = sys.argv[1:]
    keep = out = pairs = None
    name = "pr193"
    if "--keep" in args:
        i = args.index("--keep")
        if i + 1 == len(args):
            sys.exit(__doc__)
        keep = os.path.abspath(args[i + 1])
        del args[i:i + 2]
        if os.path.exists(keep) and os.listdir(keep):
            sys.exit("--keep: the directory must be empty or absent: " + keep)
    for a in [a for a in args if a.startswith("out=")]:
        out = os.path.abspath(a[4:])
        args.remove(a)
        if os.path.exists(out):
            sys.exit("out=: the file exists already: " + out)
    for a in [a for a in args if a.startswith("pairs=") or a.startswith("name=")]:
        args.remove(a)
        if a.startswith("pairs="):
            pairs = os.path.abspath(a[6:])
        else:
            name = a[5:]
    if (pairs is None) != (name == "pr193") or not re.fullmatch("[a-z0-9]+", name):
        sys.exit("pairs=<file> and name=<name> go together; <name>: lower-case letters and digits, not pr193")
    if pairs is not None and not os.path.isfile(pairs):
        sys.exit("pairs=: no such file: " + pairs)
    if len(args) != 1:
        sys.exit(__doc__)
    inputs = rbpaths.slash(args[0])
    print("step 0: sha256 of the input files in %s" % inputs, flush=True)
    if rbpaths.verify(inputs):
        sys.exit("an input file is missing or different: nothing was run.  fetch_inputs.py fetches and checks them.")
    if pairs is not None:
        got = rbpaths.sha256(pairs)
        listed = ["%s/%s" % f[:2] for f in rbpaths.OPTIONAL if f[3] == got]
        print("pairs file %s\n  sha256 %s  (%s)" % (pairs, got, "the optional input of rbpaths.py, " + listed[0]
                                                    if listed else "not a file listed in rbpaths.py"), flush=True)
    if keep is not None:
        os.makedirs(keep, exist_ok=True)
    work = rbpaths.slash(keep or tempfile.mkdtemp(prefix="rebuild-"))
    os.makedirs(work + "out")
    env = {k: w for k, w in os.environ.items() if not k.startswith("PYTHON")}
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1", REBUILD_INPUTS=inputs, REBUILD_WORK=work,
               GX_OUT=work + "gxout")
    if pairs is not None:
        env.update(REBUILD_PAIRS=pairs, REBUILD_NAME=name)
    cert = "gcert1-p11-%s.json.gz" % name
    g0, g1 = work + "out/gcert-p11-%s.json" % name, work + "out/" + cert
    steps = [["-I", "-B", HERE + "/build.py", "pr168"],
             ["-I", "-B", HERE + "/aligned.py", name],
             ["-I", "-B", HERE + "/rewrite.py", name],
             ["-I", "-B", HERE + "/phys.py", name + "xword"],
             ["-I", "-B", HERE + "/emit.py", name + "xword", g0],
             # gxconv.py finds gxpaths.py through its own directory on the module path, so it is started without -I
             ["-B", rbpaths.REPO + "tools/gx/gxconv.py", g0, "out=" + g1]]
    good = False
    try:
        for n, cmd in enumerate(steps, 1):
            print("step %d: python3 %s" % (n, " ".join(os.path.relpath(c, rbpaths.REPO) if c.endswith(".py") else
                                                         c.replace(work, "<work>/") for c in cmd)), flush=True)
            t0 = time.time()
            r = subprocess.run([sys.executable] + cmd, cwd=work, env=env)
            if r.returncode != 0:
                sys.exit("step %d failed with exit code %d" % (n, r.returncode))
            print("        (step %d: %.1f s)" % (n, time.time() - t0), flush=True)
        got, size = sha_unpacked(g1)
        print("rebuilt certificate (%d bytes, %d bytes uncompressed)" % (os.path.getsize(g1), size))
        print("  sha256 of the .gz file            %s   (differs from run to run: gzip time stamp)" % rbpaths.sha256(g1))
        print("  sha256 of the uncompressed JSON   %s" % got)
        if pairs is not None:
            print("  pairs file %s\n    sha256 %s" % (pairs, rbpaths.sha256(pairs)))
        else:
            print("  expected (tools/gx/ORIGIN.md)     %s" % EXPECTED)
        if pairs is None and os.path.isfile(rbpaths.REPO + PUBLISHED):
            print("  %s, uncompressed:" % PUBLISHED)
            print("                                    %s" % sha_unpacked(rbpaths.REPO + PUBLISHED)[0])
        if out is not None:
            os.makedirs(os.path.dirname(out), exist_ok=True)
            shutil.copyfile(g1, out)
            print("copied to " + out)
        good = got == EXPECTED or pairs is not None
    finally:
        if keep is None:
            shutil.rmtree(work)
        else:
            print("work directory kept: " + work)
    print("RESULT: %s" % ("BUILT with another pairs file, name %sxword (accepted by gx.check1 in step 6)" % name
                          if good and pairs is not None else
                          "REPRODUCED (the uncompressed content is that of the published certificate)" if good else
                          "DIFFERENT (the uncompressed content is not that of the published certificate)"))
    sys.exit(0 if good else 1)


if __name__ == "__main__":
    main()
