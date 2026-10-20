"""vcert.py: write a line-by-line zero-block certificate for named instances (rm2-solve coordinates), then CHECK it
from scratch in pure Python (no numpy): every claimed tight form is recomputed, the clique pairs are checked, and
the violating form F (Y inside D_F, |D_F| + r_F < 2|Y|) is recomputed."""
import sys, itertools, json
sys.path.insert(0, '..')
m = int(sys.argv[1]); spec = sys.argv[2]   # spec: 'ks/frames/deldir'
import os; os.environ["HM"] = str(m)
import packing2
ks = [int(v) for v in spec.split('/')[0].split(',')]; ch = [int(v) for v in spec.split('/')[1].split(',')]
z0 = int(spec.split('/')[2])
ms = [y for k, c in zip(ks, ch) for y in sorted(sorted(F) for F in packing2.all_frames(k))[c]]
ms.remove(z0)
Ts = [packing2.fd.target(k) for k in ks]
pairs = [(i, j) for i in range(m) for j in range(i, m)]
def Q(c, z): return sum(ci & ((z >> i) & 1) & ((z >> j) & 1) for ci, (i, j) in zip(c, pairs)) & 1
def p(c, T): return sum(ci & T[i][j] for ci, (i, j) in zip(c, pairs)) & 1
forms = [c for c in itertools.product([0, 1], repeat=len(pairs)) if any(c)]
info = []
for c in forms:
    D = [i for i, z in enumerate(ms) if Q(c, z)]; r = sum(p(c, T) for T in Ts)
    info.append((c, D, r))
tight = [(c, D, r) for c, D, r in info if len(D) == r and r > 0]
pairwit = {}
for c, D, r in tight:
    for i in D:
        for j in D: pairwit.setdefault((min(i, j), max(i, j)), c)
cert = None
for c, D, r in info:
    live = [i for i in D if (i, i) in pairwit]
    for sz in range(len(live), 0, -1):
        if len(D) + r >= 2 * sz: break
        for Y in itertools.combinations(live, sz):
            if all((min(a, b), max(a, b)) in pairwit for a in Y for b in Y):
                cert = dict(F=c, DF=D, rF=r, Y=list(Y), wit={f"{a},{b}": pairwit[(min(a, b), max(a, b))] for a in Y for b in Y if a <= b})
                break
        if cert: break
    if cert: break
print(f"block {m}, targets k={ks} frames {ch}, delete dir {z0}: moves (direction ints, rm2-solve coordinates) {ms}")
print("targets (rows):", [["".join(map(str, row)) for row in T] for T in Ts])
if not cert: print("NO zero-block certificate"); sys.exit()
print("CERTIFICATE (forms as coefficient vectors over monomials", pairs, ")")
print(" Y =", cert["Y"], "(moves; directions", [ms[i] for i in cert["Y"]], ")")
for key, c in cert["wit"].items(): print(f"  pair {key}: tight form {c}")
print(" F =", cert["F"], " D_F =", cert["DF"], " r_F =", cert["rF"], f" |D_F| + r_F = {len(cert['DF']) + cert['rF']} < 2|Y| = {2 * len(cert['Y'])}")
# ---- independent check, from scratch
ok = True
for key, c in cert["wit"].items():
    a, b = map(int, key.split(','))
    D = [i for i, z in enumerate(ms) if Q(c, z)]; r = sum(p(c, T) for T in Ts)
    ok &= (len(D) == r and r > 0 and a in D and b in D)
F = cert["F"]; DF = [i for i, z in enumerate(ms) if Q(F, z)]; rF = sum(p(F, T) for T in Ts)
ok &= set(cert["Y"]) <= set(DF) and len(DF) + rF < 2 * len(cert["Y"])
print("CHECK:", "PASS (the instance is excluded by the zero-block test)" if ok else "FAIL")
