#!/usr/bin/env python3
"""m-gen: assemble the Colab package of the copy-free #588 build (C588M) from the generated outputs.

pkg/c588m/overlay/          annealed variant (the build target): every NEW file, at its path in lean/complex588/
pkg/c588m/variants/{cyclic,packed}/  the files that differ for the cyclic table (fallback) and the packed table (optional upgrade)
pkg/c588m/MODULES.txt       the 70 new modules in a topological build order (imports inside the overlay), with stage
pkg/c588m/data, preflight, tools, scripts, BUILD.md, SHA256SUMS
usage: package.py <m-gen dir>"""
import sys, os, re, shutil, hashlib, json, filecmp

M = os.path.abspath(sys.argv[1])
P = os.path.join(M, "pkg", "c588m")
if os.path.exists(os.path.join(M, "pkg")):
    shutil.rmtree(os.path.join(M, "pkg"))
os.makedirs(P)


def copytree_files(src, dst):
    for root, _, files in os.walk(src):
        for f in files:
            if f.endswith(".json") and os.path.basename(root) == "logs":
                continue
            s = os.path.join(root, f)
            rel = os.path.relpath(s, src)
            if rel.startswith("logs" + os.sep):
                continue
            d = os.path.join(dst, rel)
            os.makedirs(os.path.dirname(d), exist_ok=True)
            assert not os.path.exists(d) or filecmp.cmp(s, d, shallow=False), "conflicting file %s" % rel
            shutil.copy2(s, d)


OV = os.path.join(P, "overlay")
copytree_files(os.path.join(M, "out/xmut"), OV)
copytree_files(os.path.join(M, "out/local"), OV)
copytree_files(os.path.join(M, "out/global-annealed"), OV)
os.remove(os.path.join(OV, "tools/Compare.lean"))           # identical to the tree's tools/Compare.lean: not overlaid
os.rmdir(os.path.join(OV, "tools"))
# negative comparator configuration (one unit stronger), for the official comparator if wanted
cfg = json.load(open(os.path.join(OV, "comparator/fourier-c588m.json")))
cfg["challenge_module"] = "Work.FourierC588MNeg.UniformFourierChallenge"
open(os.path.join(OV, "comparator/fourier-c588m-neg.json"), "w").write(json.dumps(cfg, indent=2) + "\n")
# variants: only the files that differ from overlay/ (cyclic = fallback, packed = optional upgrade)
nvar = {}
for var in ("cyclic", "packed"):
    VD = os.path.join(P, "variants", var)
    src = os.path.join(M, "out/global-%s" % var)
    nvar[var] = 0
    for root, _, files in os.walk(src):
        for f in files:
            s_ = os.path.join(root, f)
            rel = os.path.relpath(s_, src)
            if rel.startswith("logs") or rel.startswith("tools"):
                continue
            o = os.path.join(OV, rel)
            assert os.path.exists(o), rel
            if not filecmp.cmp(s_, o, shallow=False):
                d = os.path.join(VD, rel)
                os.makedirs(os.path.dirname(d), exist_ok=True)
                shutil.copy2(s_, d)
                nvar[var] += 1
# modules and a topological order
mods = {}
for root, _, files in os.walk(os.path.join(OV, "Work")):
    for f in files:
        if f.endswith(".lean"):
            rel = os.path.relpath(os.path.join(root, f), OV)
            mods[rel[:-5].replace(os.sep, ".")] = rel
imps = {}
for m, rel in mods.items():
    imps[m] = [x for x in re.findall(r"^import (\S+)", open(os.path.join(OV, rel), encoding="utf-8").read(), flags=re.M)
               if x in mods]
order, seen = [], set()


def visit(m, stack=()):
    assert m not in stack, "import cycle %s" % (stack,)
    if m in seen:
        return
    for x in sorted(imps[m]):
        visit(x, stack + (m,))
    seen.add(m)
    order.append(m)


def stage(m):
    if m.startswith("Work.XMut."):
        return "1-engine"
    if m.startswith("Work.Fourier") or m.endswith("Controls"):
        return "4-fourier-and-checks"
    if any(k in m for k in ("Rows", "Sched", "Inst", "Rate", "Glue401")):
        return "3-bank-rate-glue"
    return "2-certificate"


for st in ("1-engine", "2-certificate", "3-bank-rate-glue", "4-fourier-and-checks"):
    for m in sorted(mods, key=lambda x: (len(x), x)):
        if stage(m) == st:
            visit(m)
lines = ["# module  stage  (topological order: every module after all modules of the overlay that it imports)"]
lines += ["%s  %s" % (m, stage(m)) for m in order]
open(os.path.join(P, "MODULES.txt"), "w").write("\n".join(lines) + "\n")
# data, preflight, tools
D = os.path.join(P, "data")
os.makedirs(D)
for f in ["data/c588m.gcert.json.gz", "out/price/price-annealed.json", "out/price/price-annealed.schedule.json.gz",
          "out/price/price-cyclic.json", "out/price/price-cyclic.schedule.json.gz",
          "out/price/price-packed.json", "out/price/price-packed.schedule.json.gz"]:
    shutil.copy2(os.path.join(M, f), D)
PF = os.path.join(P, "preflight")
os.makedirs(PF)
for f in sorted(os.listdir(os.path.join(M, "logs"))):
    if (f.endswith(".json") and f != "package.json") or f in ("generate_local.log", "make_sched.log", "pack_leftover.log") \
            or f.startswith("generate_global"):
        shutil.copy2(os.path.join(M, "logs", f), PF)
for f in ["out/local/logs/gen-C588M.json", "out/local/logs/preflight-C588M.json"]:
    shutil.copy2(os.path.join(M, f), PF)
for v in ("annealed", "cyclic", "packed"):
    for f in ("global-C588M.json", "challenge-provenance-C588M.json"):
        shutil.copy2(os.path.join(M, "out/global-%s/logs/%s" % (v, f)), os.path.join(PF, "%s-%s" % (v, f)))
T = os.path.join(P, "tools")
shutil.copytree(os.path.join(M, "tools"), T, ignore=shutil.ignore_patterns("__pycache__"))
shutil.copy2(os.path.join(M, "REPRODUCE.sh"), T)
shutil.copytree(os.path.join(M, "scripts"), os.path.join(P, "scripts"))
shutil.copy2(os.path.join(M, "BUILD.md"), P)
# SHA256SUMS of everything
sums = []
for root, _, files in os.walk(P):
    for f in sorted(files):
        p = os.path.join(root, f)
        rel = os.path.relpath(p, P)
        if rel == "SHA256SUMS":
            continue
        sums.append("%s  %s" % (hashlib.sha256(open(p, "rb").read()).hexdigest(), rel))
sums.sort(key=lambda s: s.split("  ", 1)[1])
open(os.path.join(P, "SHA256SUMS"), "w").write("\n".join(sums) + "\n")
print(json.dumps(dict(modules=len(order), overlay_files=sum(1 for s in sums if s.split("  ")[1].startswith("overlay/")),
                      variant_files=nvar, total_files=len(sums),
                      stages={st: sum(1 for m in order if stage(m) == st) for st in
                              ("1-engine", "2-certificate", "3-bank-rate-glue", "4-fourier-and-checks")})))
