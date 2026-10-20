#!/usr/bin/env python3
"""x-words (wht45/lever2979, 19 Oct 2026): ONE stage with both #2979 words, the pair module and the generic checker,
meaning and helper files (one Lake build for everything).  Union of the two per-word stages (mk_stage.py) plus
Work/PWGen/S2979Pair.lean; refuses if the two stages disagree on a common path.  Deterministic tar (sorted, mtime 0).
usage: mk_stage_pair.py <stage dir with stage-S2979B.tgz, stage-S2979E.tgz> <pair Work dir> <out dir>"""
import sys, tarfile, io, gzip, hashlib
from pathlib import Path
sd, pairdir, od = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
files = {}
for N in ("S2979B", "S2979E"):
    with tarfile.open(sd / ("stage-%s.tgz" % N)) as tf:
        for m in tf.getmembers():
            if m.isfile():
                b = tf.extractfile(m).read(); f = m.name[2:] if m.name.startswith("./") else m.name
                assert files.get(f, b) == b, ("stages disagree", f)
                files[f] = b
for p in sorted(pairdir.rglob("*.lean")):
    f = "Work/" + str(p.relative_to(pairdir)); assert f not in files; files[f] = p.read_bytes()
buf = io.BytesIO()
with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tf:
    dirs = sorted({str(q) for f in files for q in Path(f).parents if str(q) != "."})
    for d in dirs:
        ti = tarfile.TarInfo("./" + d); ti.type = tarfile.DIRTYPE; ti.mode = 0o755; ti.mtime = 0; tf.addfile(ti)
    for f in sorted(files):
        ti = tarfile.TarInfo("./" + f); ti.size = len(files[f]); ti.mode = 0o644; ti.mtime = 0; tf.addfile(ti, io.BytesIO(files[f]))
tgz = gzip.compress(buf.getvalue(), mtime=0)
(od / "stage-S2979pair.tgz").write_bytes(tgz)
(od / "stage-S2979pair.tgz.sha256").write_text("%s  stage-S2979pair.tgz\n" % hashlib.sha256(tgz).hexdigest())
(od / "stage-S2979pair.sha256").write_text("".join("%s  ./%s\n" % (hashlib.sha256(files[f]).hexdigest(), f) for f in sorted(files)))
print("pair stage", len(files), "files; tgz", hashlib.sha256(tgz).hexdigest())
