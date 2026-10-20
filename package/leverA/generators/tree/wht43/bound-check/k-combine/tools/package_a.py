#!/usr/bin/env python3
"""k-combine: the Colab package of candidate (A), k-lower's U8X (exit frames only, 1.2704023e-3), re-assembled with
the x-exit-lean modules, the negative comparator config, generic scripts and BUILD.md.  The overlay files are
k-lower's package files unchanged (sha256 compared with k-lower's SHA256SUMS.overlay).
usage: package_a.py <k-combine dir> <k-lower pkg dir> <x-exit-lean XExit dir> <k-lower colab scripts dir>"""
import sys, os, hashlib, shutil, json, subprocess
from pathlib import Path
K, KP, XE, KS = map(Path, sys.argv[1:5])
NAME = "U8X"
pkg = K / "pkg" / "k-combine-u8x"
if pkg.exists(): shutil.rmtree(pkg)
shutil.copytree(KP / "overlay", pkg / "overlay")
# k-lower's checksums still hold
for line in (KP / "SHA256SUMS.overlay").read_text().split("\n"):
    if not line.strip(): continue
    h, rel = line.split("  ", 1)
    assert hashlib.sha256((pkg / rel).read_bytes()).hexdigest() == h, rel
for f in sorted(XE.glob("*.lean")):
    (pkg / "overlay/XExit").mkdir(parents=True, exist_ok=True); shutil.copy2(f, pkg / "overlay/XExit" / f.name)
c = pkg / "overlay/comparator"
pos = json.loads((c / "fourier-u8x.json").read_text())
(c / "fourier-u8x-neg.json").write_text(json.dumps(dict(pos, challenge_module="Work.FourierU8XNeg.UniformFourierChallenge"), indent=2) + "\n")
shutil.copy2(KP / "MODULES.txt", pkg / "MODULES.txt")
lines = ["%s  overlay/%s" % (hashlib.sha256(p.read_bytes()).hexdigest(), p.relative_to(pkg / "overlay"))
         for p in sorted((pkg / "overlay").rglob("*")) if p.is_file()]
(pkg / "SHA256SUMS.overlay").write_text("\n".join(lines) + "\n")
sc = pkg / "scripts"; sc.mkdir()
for t in ("setup.sh", "build-direct.sh", "cmp.sh"):
    shutil.copy2(K / "tools/colab" / t, sc / t); os.chmod(sc / t, 0o755)
shutil.copy2(KS / "checks.sh", sc / "checks.sh"); os.chmod(sc / "checks.sh", 0o755)
if (K / "docs/BUILD-U8X.md").exists(): shutil.copy2(K / "docs/BUILD-U8X.md", pkg / "BUILD.md")
allsum = ["%s  %s" % (hashlib.sha256(p.read_bytes()).hexdigest(), p.relative_to(pkg))
          for p in sorted(pkg.rglob("*")) if p.is_file() and p.name != "SHA256SUMS"]
(pkg / "SHA256SUMS").write_text("\n".join(allsum) + "\n")
tb = str(pkg) + ".tar.gz"
subprocess.run(["tar", "--no-xattrs", "-czf", tb, "-C", str(pkg.parent), pkg.name], check=True, env=dict(os.environ, COPYFILE_DISABLE="1"))
h = hashlib.sha256(open(tb, "rb").read()).hexdigest()
open(tb + ".sha256", "w").write("%s  %s\n" % (h, os.path.basename(tb)))
print(json.dumps(dict(files=len(lines), modules=len((pkg / "MODULES.txt").read_text().split()), tarball=tb, sha256=h)))
