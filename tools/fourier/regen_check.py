#!/usr/bin/env python3
"""Regenerate the Fourier copies from the unmodified upstream files and compare bytes.

usage (from the repository root):  python3 tools/fourier/regen_check.py [<seam module> = Work.Fourier.Seam]

tools/fourier/mkchain.py writes, into a temporary directory, the eleven modified copies of the files of
OAI/Computability/FourierTransform/, the file Goal.lean, the comparator challenge and its configuration.  Each
generated file is compared byte for byte with the file of the same name in Work/Fourier/ (the configuration:
with every file in comparator/).  A generated file whose name is not in the repository is compared with every
file of Work/Fourier/, so that a renamed file is still found.  Exit status 0 only if nothing differs.
"""
import filecmp
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
seam = sys.argv[1] if len(sys.argv) > 1 else "Work.Fourier.Seam"
out = tempfile.mkdtemp(prefix="fourier-regen-")
try:
    r = subprocess.run([sys.executable, "-B", os.path.join(ROOT, "tools/fourier/mkchain.py"), seam],
                       env=dict(os.environ, FOURIER_OUT=out), capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode != 0:
        print(r.stderr.strip())
        sys.exit("RESULT: FAILED (mkchain.py exit %d)" % r.returncode)
    same = diff = absent = 0
    repo = os.path.join(ROOT, "Work/Fourier")
    for f in sorted(os.listdir(out)):
        g = os.path.join(out, f)
        if f.endswith(".json"):
            if f == "rename-map.json":
                continue
            pool = [os.path.join(ROOT, "comparator", x) for x in sorted(os.listdir(os.path.join(ROOT, "comparator")))]
            hit = [p for p in pool if filecmp.cmp(g, p, shallow=False)]
            print("%-40s %s" % (f, "identical to comparator/" + os.path.basename(hit[0]) if hit else "(no identical file in comparator/)"))
            same, absent = same + bool(hit), absent + (not hit)
            continue
        p = os.path.join(repo, f)
        if os.path.isfile(p):
            ok = filecmp.cmp(g, p, shallow=False)
            print("%-40s %s" % (f, "identical" if ok else "DIFFERENT"))
            same, diff = same + ok, diff + (not ok)
        else:
            hit = [x for x in sorted(os.listdir(repo)) if filecmp.cmp(g, os.path.join(repo, x), shallow=False)]
            print("%-40s %s" % (f, "identical to Work/Fourier/" + hit[0] if hit else "(not in the repository)"))
            same, absent = same + bool(hit), absent + (not hit)
    print("RESULT: %s (%d repository files reproduced byte for byte, %d different, %d generated files not in the repository)"
          % ("REPRODUCED" if diff == 0 and same > 0 else "NOT REPRODUCED", same, diff, absent))
    sys.exit(0 if diff == 0 and same > 0 else 1)
finally:
    shutil.rmtree(out, ignore_errors=True)
