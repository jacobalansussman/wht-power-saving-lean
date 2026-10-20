import itertools, numpy as np
from vlib import *
exec(open("v2_f21.py").read().split("# all elements g")[0].split("c2 = c @ c")[1])  # dq, SYMS, fixscalar
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
c = np.load("v_c.npy"); n = np.load("v_n.npy"); EN = np.load("v_EN.npy")
L = lagrangians(); d1 = [j for j in range(len(L)) if len(L[j] & LZ) == 4]
sn = sym(n); fix1 = [L[j] for j in d1 if apply_sym(sn, L[j]) == L[j]]
say("n-fixed distance-1 cells:", [sorted(x) for x in fix1])
reps = {}
for Sm in SYMS:
    D = dq(Sm)
    for A in gl3():
        PA = permm(lambda x, A=A: lin(A, x))
        for z in range(1, 8):
            m = D @ PA @ dirm(z)
            cl = apply_sym(sym(m), LZ)
            if cl in fix1 and cl not in reps: reps[cl] = m
        if len(reps) == 2: break
    if len(reps) == 2: break
say("reps found:", len(reps))
def sigma_of(gm, m):
    NM = orth(np.hstack([Nb, rho(m) @ Nb]))
    Rg = [rho(np.linalg.matrix_power(gm, s)) for s in range(3)]
    assert np.linalg.matrix_rank(np.hstack([NM, Rg[1] @ NM]), tol=1e-7) == NM.shape[1]
    return mults([quot_trace(R, NM, Nb) for R in Rg], 3)
tw = []
for a in range(8):
    for b in range(8):
        g = fixscalar(n @ pau(a, b), 3)
        if g is not None: tw.append(((a, b), g))
for cl, m in reps.items():
    say("cell", sorted(cl), " sigma(n) =", sigma_of(n, m))
res = {}
for ab, g in tw:
    dd = mults([quot_trace(rho(np.linalg.matrix_power(Kinv @ inv(g) @ Kmat, s)), Nb, EN) for s in range(3)], 3)
    ss = tuple(tuple(sigma_of(g, m)) for cl, m in reps.items())
    res.setdefault((tuple(dd), ss), []).append(ab)
for k, v in res.items(): say("lift n X^a Z^b: demand split", k[0], " cap splits on the two cells", k[1], " count", len(v), v[:4])
# exact commuters in the second cell, chunked over B
DS = np.array([dq(Sm) @ pau(a, b) for Sm in SYMS for a in range(8) for b in range(8)])
for cl, m in reps.items():
    tot = [0, 0]
    for A in gl3():
        PA = permm(lambda x, A=A: lin(A, x))
        mm = m[None] @ (PA[None] @ DS)
        fl = (n[None] @ mm).reshape(len(mm), -1); fr = (mm @ n[None]).reshape(len(mm), -1)
        k0 = np.argmax(np.abs(fr), axis=1)
        ratio = fl[np.arange(len(mm)), k0] / fr[np.arange(len(mm)), k0]
        ok = np.all(np.abs(fl - ratio[:, None] * fr) < 1e-8, axis=1)
        tot[0] += int(ok.sum()); tot[1] += int((ok & (np.abs(ratio - 1) < 1e-8)).sum())
    say("cell", sorted(cl), ": letters (mod scalar) commuting with n up to scalar", tot[0], " exactly", tot[1])
open("out_v3_cells.txt", "w").write("\n".join(out) + "\n")
