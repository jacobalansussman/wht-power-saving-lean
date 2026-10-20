#!/usr/bin/env python3
"""p-wgen: the stage of one warped word overlay: the generated Work/ tree of out/gen-<NAME> (data, checks, price, rate,
schedule, glue, Fourier files, challenges), the comparator configurations (positive, and negative with the Neg
challenge module), p-wgen's Work/PWGen/Fp.lean and the checker files it was built against (Work/KLowW/*.lean, copied
from p-warpchk; recorded by sha256).  Writes stage-<NAME>.tgz (+ .tgz.sha256), stage-<NAME>.sha256 (per file) and
stage-<NAME>-modules.txt (all modules of the stage).  Deterministic tar (sorted, mtime 0, uid 0).
usage: mk_stage.py <NAME> <gen dir> <extra Work dir (KLowW, Warp, PWGen)> <unused> <out dir>"""
import sys, os, json, tarfile, hashlib, io, gzip
from pathlib import Path

NAME, gen, klw, fp, od = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]), Path(sys.argv[5])
files = {}
for p in sorted((gen / "Work").rglob("*.lean")):
    files["Work/" + str(p.relative_to(gen / "Work"))] = p.read_bytes()
for p in sorted(klw.rglob("*.lean")):             # the checker/meaning/helper files the stage was built with
    files["Work/" + str(p.relative_to(klw))] = p.read_bytes()
cp = gen / "comparator" / ("fourier-%s.json" % NAME.lower())
pos = json.loads(cp.read_text())
neg = dict(pos, challenge_module="Work.Fourier%sNeg.UniformFourierChallenge" % NAME)
files["comparator/fourier-%s.json" % NAME.lower()] = cp.read_bytes()
files["comparator/fourier-%s-neg.json" % NAME.lower()] = (json.dumps(neg, indent=2) + "\n").encode()
od.mkdir(parents=True, exist_ok=True)
buf = io.BytesIO()
with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tf:
    dirs = sorted({str(Path(f).parent) for f in files} | {str(q) for f in files for q in Path(f).parents if str(q) != "."})
    for d in dirs:
        ti = tarfile.TarInfo("./" + d); ti.type = tarfile.DIRTYPE; ti.mode = 0o755; ti.mtime = 0
        tf.addfile(ti)
    for f in sorted(files):
        ti = tarfile.TarInfo("./" + f); ti.size = len(files[f]); ti.mode = 0o644; ti.mtime = 0
        tf.addfile(ti, io.BytesIO(files[f]))
tgz = gzip.compress(buf.getvalue(), mtime=0)
(od / ("stage-%s.tgz" % NAME)).write_bytes(tgz)
(od / ("stage-%s.tgz.sha256" % NAME)).write_text("%s  stage-%s.tgz\n" % (hashlib.sha256(tgz).hexdigest(), NAME))
(od / ("stage-%s.sha256" % NAME)).write_text("".join("%s  ./%s\n" % (hashlib.sha256(files[f]).hexdigest(), f) for f in sorted(files)))
mods = sorted(f[:-5].replace("/", ".") for f in files if f.endswith(".lean"))
(od / ("stage-%s-modules.txt" % NAME)).write_text("\n".join(mods) + "\n")
ax = ["import Work.GCert.Data.%sGlue401" % NAME, "import Work.PWGen.B401W.EndW"] + ["#print axioms " + x for x in [
    "GXD.%s.ginvW" % NAME, "GXD.%s.ginvXW" % NAME, "GXD.%s.ginv_cost" % NAME, "GXD.%s.price" % NAME,
    "GXD.%s.fitsS" % NAME, "GXD.%s.hent" % NAME, "GXD.%s.width" % NAME, "GXD.%s.okC" % NAME, "GXD.%s.wvalid" % NAME,
    "GXD.%s.res_ok" % NAME, "GXD.%s.fp_ok" % NAME, "GXD.%s.rank1_ok" % NAME, "GXD.%s.sched_ok" % NAME,
    "OAI.PowerSaving.BridgeRate.%s.unit_cost_x" % NAME, "OAI.PowerSaving.BlockAccounting.fold_%sx" % NAME,
    "PWG.roles", "PWG.res", "OAI.PowerSaving.GX.engine_of_ginv_401W", "KLW.WValid.invW"]] + \
    ["#check @GXD.%s.ginvW" % NAME, "#check @GXD.%s.ginv_cost" % NAME]
(od / ("Ax%s.lean" % NAME)).write_text("\n".join(ax) + "\n")
print(NAME, len(files), "files", sum(len(b) for b in files.values()), "bytes; tgz", hashlib.sha256(tgz).hexdigest())
