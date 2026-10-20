#!/usr/bin/env python3
"""k-lower: assemble the Colab overlay (new files only, paths relative to lean/complex588), the topological module
list, SHA256SUMS.  usage: package.py <k-lower dir> <NAME>"""
import sys, os, re, hashlib, shutil, json
from pathlib import Path
import argparse
ap = argparse.ArgumentParser(); ap.add_argument("k"); ap.add_argument("name")
ap.add_argument("--local", default="out/local"); ap.add_argument("--global", dest="glob", default="out/global")
a = ap.parse_args()
K = Path(a.k); NAME = a.name
pkg = K / "pkg" / ("k-lower-%s" % NAME.lower())
if pkg.exists(): shutil.rmtree(pkg)
ov = pkg / "overlay"
src = []
for f in sorted((K / "proj/Work/KLow").glob("*.lean")): src.append((f, Path("Work/KLow") / f.name))
for root in (K / a.local, K / a.glob):
    for f in sorted(root.rglob("*.lean")):
        rel = f.relative_to(root)
        if rel.parts[0] == "Work": src.append((f, rel))
for f in sorted((K / a.glob / "comparator").glob("*.json")): src.append((f, Path("comparator") / f.name))
mods = {}
for f, rel in src:
    d = ov / rel; d.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(f, d)
    if rel.suffix == ".lean":
        mods[".".join(rel.with_suffix("").parts)] = d
deps = {m: [i for i in re.findall(r"^import (\S+)", p.read_text(), re.M) if i in mods] for m, p in mods.items()}
order, seen = [], set()
def visit(m):
    if m in seen: return
    seen.add(m)
    for d in deps[m]: visit(d)
    order.append(m)
for m in sorted(mods): visit(m)
skip = {"Work.Fourier%s.UniformFourierChallenge" % NAME, "Work.Fourier%sNeg.UniformFourierChallenge" % NAME, "Work.Fourier%s.Axioms" % NAME}
(pkg / "MODULES.txt").write_text("\n".join(m for m in order if m not in skip) + "\n")
lines = []
for p in sorted(ov.rglob("*")):
    if p.is_file():
        lines.append("%s  %s" % (hashlib.sha256(p.read_bytes()).hexdigest(), p.relative_to(ov)))
(pkg / "SHA256SUMS.overlay").write_text("\n".join("%s  overlay/%s" % tuple(l.split("  ", 1)) for l in lines) + "\n")
sc = pkg / "scripts"; sc.mkdir()
for t in ("setup.sh", "build.sh"):
    txt = (K / "tools/colab" / t).read_text().replace("@NAME@", NAME).replace("@name@", NAME.lower())
    (sc / t).write_text(txt); os.chmod(sc / t, 0o755)
import subprocess
subprocess.run(["tar", "--no-xattrs", "-czf", str(pkg) + ".tar.gz", "-C", str(pkg.parent), pkg.name], check=True,
               env=dict(os.environ, COPYFILE_DISABLE="1"))
print(len(mods), "modules,", len(lines), "files ->", pkg)
