#!/usr/bin/env python3
"""p-wgen: mutants of a gcert/1w word; tools/wgen2.py --mirror-only must REJECT each one.
usage: mutants_w.py <word.gcert1w.json.gz> <scratch dir>"""
import sys, json, gzip, copy, subprocess, os
c = json.load(gzip.open(sys.argv[1])); od = sys.argv[2]; EXP = sys.argv[3]     # x-words: expected profile
here = os.path.dirname(os.path.abspath(__file__))
def m_gcol(d):      # a wrong G column (still invertible: swap two columns)
    g = d["warps"][0]["G_columns"]; g[0], g[1] = g[1], g[0]
def m_unwarp(d):    # drop one unwarp of U
    i = next(i for i, g in enumerate(d["U"]) if g[0] == "warp"); del d["U"][i]
def m_coef(d):      # double one compensation coefficient
    g = next(g for g in d["C0"] if g[0] == "out" and g[3]); g[3][0][1] *= 2
def m_tadd(d):      # flip the sign of one twisted-add warp letter
    g = next(g for g in d["D0"] if g[0] == "warp"); g[2] = -g[2]
def m_frame(d):     # a B warp ends in the wrong frame
    g = next(g for g in d["B"] if g[0] == "warp"); g[4] = g[3]
def m_scat(d):      # change one scatter coefficient
    row = next(r for r in d["scat"]["table"] if r); row[0][1] += 1
def m_zero(d):      # drop the last READ of a zero buffer into a live register (z -> x/slot): the live map changes
    i = max(i for i, g in enumerate(d["U"]) if g[0] == "out" and g[2] >= 632 and any(t < 632 for t, _, _ in g[3])); del d["U"][i]
# NOT a mutant in our model: dropping the last WRITE into a zero buffer (the hub's "zero at zfree") leaves the live
# map unchanged; `.zf` erases the private scratch (p-warpchk's route semantics), so the generator accepts it.
res = {}
for nm, f in [("wrong G column", m_gcol), ("unwarp dropped", m_unwarp), ("compensation coefficient doubled", m_coef),
              ("twisted-add letter flipped", m_tadd), ("B warp end frame wrong", m_frame), ("scatter coefficient changed", m_scat),
              ("zero-buffer read dropped", m_zero)]:
    d = copy.deepcopy(c); f(d)
    p = os.path.join(od, "mut.json.gz"); open(p, "wb").write(gzip.compress(json.dumps(d).encode(), mtime=0))
    r = subprocess.run([sys.executable, "-I", os.path.join(here, "wgen2.py"), p, "MUT", od, "--mirror-only", "--expect", EXP], capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip().splitlines()
    res[nm] = (r.returncode, out[-1][:160] if out else "")
    print(nm, "->", "REJECTED" if r.returncode == 2 else "ACCEPTED (BAD)", res[nm][1])
assert all(rc == 2 for rc, _ in res.values()), "a mutant was accepted"
print("all mutants rejected")
