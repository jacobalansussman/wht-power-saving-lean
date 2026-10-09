#!/usr/bin/env python3
"""Regenerate the generated Lean files from the certificate and compare them with the files of this repository.

usage (from anywhere):  python3 tools/regen_check.py [--keep <dir>]

Runs, into a fresh directory (a temporary one, or <dir> with --keep; never into the source tree):
    tools/gen/ccgen.py tools/certificate/c2-combine-best-h16.json Cr2h16
    tools/gen/crate.py B2Ke16 Cr2h16
and compares every generated file, byte for byte, with the file of the same path in this repository.
Exit code 0 = the certificate has the expected sha256 and all 16 files that are part of the repository are
reproduced byte for byte.  The generators also write five files that are NOT part of the repository (the
eight-decimal comparator package and two information logs); they are listed, not compared.

Needs only Python 3 (standard library).  No Lean is run.
"""
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CERT = "tools/certificate/c2-combine-best-h16.json"
CERT_SHA256 = "1325bd8a6ccf1bfa2b773d9846c2e7dc8c120c4860c2b4317e109de3e3a20d51"
CERT_NAME, INST = "Cr2h16", "B2Ke16"

IN_REPO = (["Work/CarrierCheck/Gen/%s%s.lean" % (CERT_NAME, s) for s in ("", "Lab", "Shape", "Scal", "Hist")]
           + ["Work/CarrierCheck/%sCert.lean" % CERT_NAME]
           + ["Work/CarrierCheck/%s%s.lean" % (p, INST + s) for p, s in
              (("Rate", ""), ("Rate", "x"), ("Inst", ""), ("", ""),
               ("Challenge", "x"), ("Solution", "x"), ("Challenge", "xR"), ("Solution", "xR"))]
           + ["comparator/%sx.json" % INST, "comparator/%sxR.json" % INST])
NOT_IN_REPO = ["Work/CarrierCheck/Challenge%s.lean" % INST, "Work/CarrierCheck/Solution%s.lean" % INST,
               "comparator/%s.json" % INST, "tools/gen/logs/rate-%s.json" % INST, "tools/gen/logs/rate-%sx.json" % INST]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    keep = None
    if "--keep" in sys.argv:
        keep = os.path.abspath(sys.argv[sys.argv.index("--keep") + 1])
        if os.path.exists(keep) and os.listdir(keep):
            sys.exit("--keep: the directory must be empty or absent: " + keep)
        os.makedirs(keep, exist_ok=True)
    out = keep or tempfile.mkdtemp(prefix="regen-")
    bad = same = 0
    got = sha(os.path.join(ROOT, CERT))
    print("certificate %s\n  sha256 %s  %s" % (CERT, got, "as expected" if got == CERT_SHA256 else "NOT THE EXPECTED ONE"))
    bad += got != CERT_SHA256
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    for cmd in (["tools/gen/ccgen.py", CERT, CERT_NAME], ["tools/gen/crate.py", INST, CERT_NAME]):
        print("$ python3 " + " ".join(cmd) + " --out <fresh directory>", flush=True)
        r = subprocess.run([sys.executable] + cmd + ["--out", out], cwd=ROOT, env=env)
        if r.returncode != 0:
            sys.exit("generator failed with exit code %d" % r.returncode)
    made = sorted(os.path.relpath(os.path.join(d, f), out) for d, _, fs in os.walk(out) for f in fs)
    for rel in IN_REPO:
        a, b = os.path.join(out, rel), os.path.join(ROOT, rel)
        if not os.path.isfile(a):
            res = "NOT GENERATED"
        elif not os.path.isfile(b):
            res = "MISSING IN THE REPOSITORY"
        else:
            res = "IDENTICAL" if open(a, "rb").read() == open(b, "rb").read() else "DIFFERENT"
        bad += res != "IDENTICAL"
        same += res == "IDENTICAL"
        print("%-9s %9d bytes  %s  %s" % (res, os.path.getsize(a) if os.path.isfile(a) else 0,
                                           sha(a)[:16] if os.path.isfile(a) else "-" * 16, rel))
    for rel in NOT_IN_REPO:
        print("%-9s %s" % ("(extra)" if rel in made else "NOT GENERATED", rel + "   -- generated, not part of this repository"))
    unexpected = [m for m in made if m not in IN_REPO and m not in NOT_IN_REPO]
    for m in unexpected:
        print("UNEXPECTED OUTPUT " + m)
    bad += len(unexpected)
    if keep is None:
        shutil.rmtree(out)
    else:
        print("generated files kept in " + keep)
    print("RESULT: %s (%d of %d repository files reproduced byte for byte)"
          % ("REPRODUCED" if bad == 0 else "FAIL", same, len(IN_REPO)))
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
