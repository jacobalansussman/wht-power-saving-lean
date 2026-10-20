#!/usr/bin/env python3
"""k-combine: assemble the Colab overlay of the stacked certificate NAME (new files only, paths relative to
lean/complex588): x-mutable's engine Work/XMut (m-gen's relocation, byte-identical to the C588M build), k-lower's
Work/KLow, k-combine's Work/KMX (Colab variant), the generated local and global modules, the comparator configs
(positive and one-unit-stronger negative); MODULES.txt in topological order (challenges, Axioms and Controls are
built by checks.sh); SHA256SUMS of the overlay; the scripts; the tarball and its own sha256.
usage: package_xz.py <k-combine dir> <NAME> <m-gen XMut overlay dir> <k-lower KLow overlay dir>"""
import sys, os, re, hashlib, shutil, json, subprocess
from pathlib import Path
K, NAME, XM, KL = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
pkg = K / "pkg" / ("k-combine-%s" % NAME.lower())
if pkg.exists(): shutil.rmtree(pkg)
ov = pkg / "overlay"
src = []
for f in sorted(XM.glob("*.lean")): src.append((f, Path("Work/XMut") / f.name))
for f in sorted(KL.glob("*.lean")): src.append((f, Path("Work/KLow") / f.name))
for f in sorted((K / "src-colab/Work/KMX").glob("*.lean")): src.append((f, Path("Work/KMX") / f.name))
XE = K.resolve().parent / "x-exit-lean/proj/XExit"
assert len(list(XE.glob("*.lean"))) == 5, XE
for f in sorted(XE.glob("*.lean")): src.append((f, Path("XExit") / f.name))
for root in (K / ("out/local-%s" % NAME), K / ("out/global-%s" % NAME)):
    for f in sorted(root.rglob("*.lean")):
        rel = f.relative_to(root)
        if rel.parts[0] == "Work": src.append((f, rel))
g = K / ("out/global-%s" % NAME) / "comparator"
pos = json.loads((g / ("fourier-%s.json" % NAME.lower())).read_text())
neg = dict(pos, challenge_module="Work.Fourier%sNeg.UniformFourierChallenge" % NAME)
(g / ("fourier-%s-neg.json" % NAME.lower())).write_text(json.dumps(neg, indent=2) + "\n")
for f in sorted(g.glob("*.json")): src.append((f, Path("comparator") / f.name))
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
skip = {"Work.Fourier%s.UniformFourierChallenge" % NAME, "Work.Fourier%sNeg.UniformFourierChallenge" % NAME,
        "Work.Fourier%s.Axioms" % NAME, "Work.GCert.Data.%sControls" % NAME} | {m for m in mods if m.startswith("XExit.")}
(pkg / "MODULES.txt").write_text("\n".join(m for m in order if m not in skip) + "\n")
lines = []
for p in sorted(ov.rglob("*")):
    if p.is_file():
        lines.append("%s  overlay/%s" % (hashlib.sha256(p.read_bytes()).hexdigest(), p.relative_to(ov)))
(pkg / "SHA256SUMS.overlay").write_text("\n".join(lines) + "\n")
sc = pkg / "scripts"; sc.mkdir()
for t in ("setup.sh", "build-direct.sh", "checks.sh", "cmp.sh"):
    shutil.copy2(K / "tools/colab" / t, sc / t); os.chmod(sc / t, 0o755)
for extra in ("BUILD.md",):
    if (K / "docs" / ("BUILD-%s.md" % NAME)).exists(): shutil.copy2(K / "docs" / ("BUILD-%s.md" % NAME), pkg / extra)
allsum = []
for p in sorted(pkg.rglob("*")):
    if p.is_file() and p.name != "SHA256SUMS":
        allsum.append("%s  %s" % (hashlib.sha256(p.read_bytes()).hexdigest(), p.relative_to(pkg)))
(pkg / "SHA256SUMS").write_text("\n".join(allsum) + "\n")
tb = str(pkg) + ".tar.gz"
subprocess.run(["tar", "--no-xattrs", "-czf", tb, "-C", str(pkg.parent), pkg.name], check=True,
               env=dict(os.environ, COPYFILE_DISABLE="1"))
h = hashlib.sha256(open(tb, "rb").read()).hexdigest()
open(tb + ".sha256", "w").write("%s  %s\n" % (h, os.path.basename(tb)))
print(json.dumps(dict(modules=len(mods), built_by_build_direct=len([m for m in order if m not in skip]), files=len(lines),
                      package=str(pkg), tarball=tb, sha256=h)))
