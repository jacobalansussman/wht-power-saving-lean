"""term-rank certificate for the rm-theory set's 8 in 39 (rm2-solve chain hR1: plain 8 in 40, order rand1, delete
move 1 = direction 4), in rm2-solve's coordinates. Finds the form F with term rank < |D_F| - r_F, and a Konig cover
(rows Rc, columns Cc covering every entry not forced to zero), then re-checks everything from scratch."""
import sys, itertools
import numpy as np
sys.path.insert(0, '..')
from inst import canonical
ks = [0, 2, 4, 8, 18, 24, 25, 27]
inst, fd = canonical(5, ks, dirs="frames", order="rand1", seed=0)
Ts = [np.array(T) for T in inst.T]; dirs = list(inst.dirs)
assert dirs[1] == 4
ms = dirs[:1] + dirs[2:]
m = 5; pairs = [(i, j) for i in range(m) for j in range(i, m)]
C = np.array(list(itertools.product([0, 1], repeat=len(pairs))))[1:]
Z = np.array([[(z >> j) & 1 for j in range(m)] for z in ms])
mon = np.stack([Z[:, i] * Z[:, j] for (i, j) in pairs], 1)
Emv = (C @ mon.T) % 2
tv = np.stack([[T[i, j] for (i, j) in pairs] for T in Ts], 1)
r = ((C @ tv) % 2).sum(1); D = Emv.sum(1)
assert (D >= r).all()
tight = (D == r) & (r > 0)
Tm = Emv[tight]; fz = (Tm.T @ Tm) > 0
print("8 in 39 moves:", ms); print("tight forms:", int(tight.sum()), " live moves:", [i for i in range(len(ms)) if fz[i, i]])
def matching(Dset):
    allowed = {i: [j for j in Dset if not fz[i, j]] for i in Dset}; match = {}
    def aug(i, vis):
        for j in allowed[i]:
            if j in vis: continue
            vis.add(j)
            if j not in match or aug(match[j], vis): match[j] = i; return True
        return False
    for i in Dset: aug(i, set())
    return match, allowed
found = None
for f in range(len(C)):
    Dset = list(np.nonzero(Emv[f])[0])
    if not fz[np.ix_(Dset, Dset)].any(): continue
    match, allowed = matching(Dset)
    if len(match) < D[f] - r[f]: found = (f, Dset, match, allowed); break
f, Dset, match, allowed = found
# Konig cover from the matching: Z0 = rows unmatched + alternating reachable; cover = (rows not in Z) + (cols in Z)
mrow = {i: j for j, i in match.items()}
visR, visC = set(), set(); stack = [i for i in Dset if i not in mrow]
while stack:
    i = stack.pop()
    if i in visR: continue
    visR.add(i)
    for j in allowed[i]:
        if j not in visC:
            visC.add(j)
            if j in match: stack.append(match[j])
Rc = [i for i in Dset if i not in visR]; Cc = sorted(visC)
print("F =", tuple(C[f]), " D_F =", Dset, f"(directions {[ms[i] for i in Dset]})", " r_F =", int(r[f]))
print(f"need rank E[D_F,D_F] = |D_F| - r_F = {len(Dset) - int(r[f])}; term rank = {len(match)}")
print("Konig cover: rows", Rc, " columns", Cc, f" size {len(Rc) + len(Cc)}")
need = {(i, j) for i in Dset for j in Dset if fz[i, j]}
wit = {}
for (i, j) in need:
    for g in np.nonzero(tight)[0]:
        if Emv[g, i] and Emv[g, j]: wit[(i, j)] = tuple(C[g]); break
print("forced zeros inside D_F (pairs sharing a tight support):", sorted(need))
for (i, j) in sorted(need):
    if i <= j: print(f"  ({i},{j}) witness tight form {wit[(i, j)]}")
# --- from-scratch re-check in pure Python
def Qv(c, z): return sum(ci & ((z >> a) & 1) & ((z >> b) & 1) for ci, (a, b) in zip(c, pairs)) & 1
def pv(c, T): return sum(ci & int(T[a][b]) for ci, (a, b) in zip(c, pairs)) & 1
ok = True
for (i, j), c in wit.items():
    Dg = [t for t, z in enumerate(ms) if Qv(c, z)]; rg = sum(pv(c, T) for T in Ts)
    ok &= len(Dg) == rg > 0 and i in Dg and j in Dg
cF = tuple(C[f]); DF = [t for t, z in enumerate(ms) if Qv(cF, z)]; rF = sum(pv(cF, T) for T in Ts)
ok &= DF == Dset and rF == int(r[f])
cover_ok = all((i in Rc) or (j in Cc) for i in DF for j in DF if (i, j) not in need)
ok &= cover_ok and len(Rc) + len(Cc) < len(DF) - rF
print("CHECK:", "PASS: rank E[D_F,D_F] <= cover size < |D_F| - r_F, so this 8 in 39 is impossible" if ok else "FAIL")
