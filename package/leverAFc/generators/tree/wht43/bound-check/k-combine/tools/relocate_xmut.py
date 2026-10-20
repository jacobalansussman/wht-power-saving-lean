#!/usr/bin/env python3
"""m-gen: x-mutable's 7 engine modules (x-mutable/proj/XMut/*.lean, compiled on the Mac 13 Oct, only the 3 standard
axioms) placed at Work/XMut/*.lean for Shea's tree, whose lakefile builds `lean_lib Work` (globs Work.+).
The ONLY change: `import XMut.<M>` -> `import Work.XMut.<M>` (namespaces are untouched: XMut, GS, GX, GLab, GXD).
Writes <out>/Work/XMut/*.lean and prints the sha256 of every source and output.
usage: relocate_xmut.py <x-mutable/proj/XMut> <out>"""
import sys, os, re, hashlib, json
src, out = sys.argv[1], sys.argv[2]
os.makedirs(os.path.join(out, "Work/XMut"), exist_ok=True)
rec = {}
for f in sorted(os.listdir(src)):
    if not f.endswith(".lean"): continue
    t = open(os.path.join(src, f), encoding="utf-8").read()
    n = len(re.findall(r"^import XMut\.", t, flags=re.M))
    u = re.sub(r"^import XMut\.", "import Work.XMut.", t, flags=re.M)
    assert "XMut.lean" not in u and not re.search(r"^import XMut", u, flags=re.M)
    open(os.path.join(out, "Work/XMut", f), "w", encoding="utf-8").write(u)
    rec[f] = dict(source_sha256=hashlib.sha256(t.encode()).hexdigest(), out_sha256=hashlib.sha256(u.encode()).hexdigest(),
                  imports_rewritten=n, lines=t.count("\n"))
print(json.dumps(rec, indent=1))
