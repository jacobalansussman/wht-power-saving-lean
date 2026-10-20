import itertools, numpy as np
from vlib import *
exec(open("v2_f21.py").read().split("# all elements g")[0].split("c2 = c @ c")[1])  # dq, SYMS, fixscalar
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
L = lagrangians(); d1 = [j for j in range(len(L)) if len(L[j] & LZ) == 4]
def bitperm(p):  # coordinate permutation: bit j of x goes to bit p[j]
    def f(x): return sum(((x >> j) & 1) << p[j] for j in range(3))
    return permm(f)
t3 = bitperm((1, 2, 0)); t2 = bitperm((1, 0, 2))
for nm, g in [("3-cycle", t3), ("transposition (01)", t2)]:
    assert np.allclose(g @ Kmat, Kmat @ g)
    sg = sym(g); fx = [L[j] for j in d1 if apply_sym(sg, L[j]) == L[j]]
    say(nm, ": fixed distance-1 cells:", len(fx))
DS = np.array([dq(Sm) @ pau(a, b) for Sm in SYMS for a in range(8) for b in range(8)])
sg3 = sym(t3); fx3 = [L[j] for j in d1 if apply_sym(sg3, L[j]) == L[j]]
reps = {}
for Sm in SYMS:
    D = dq(Sm)
    for A in gl3():
        PA = permm(lambda x, A=A: lin(A, x))
        for z in range(1, 8):
            m = D @ PA @ dirm(z); cl = apply_sym(sym(m), LZ)
            if cl in fx3 and cl not in reps: reps[cl] = m
        if len(reps) == len(fx3): break
    if len(reps) == len(fx3): break
say("reps for 3-cycle-fixed cells:", len(reps))
def commute_count(g, m, also=None):
    tot = [0, 0]
    for A in gl3():
        PA = permm(lambda x, A=A: lin(A, x))
        mm = m[None] @ (PA[None] @ DS)
        ok = np.all(np.abs((g[None] @ mm) - (mm @ g[None])) < 1e-8, axis=(1, 2))
        tot[0] += int(ok.sum())
        if also is not None:
            ok2 = ok & np.all(np.abs((also[None] @ mm) - (mm @ also[None])) < 1e-8, axis=(1, 2))
            tot[1] += int(ok2.sum())
    return tot
for cl, m in reps.items():
    a3, both = commute_count(t3, m, also=t2)
    say("3-cycle-fixed cell", sorted(cl), ": letters (mod scalar) commuting exactly with the 3-cycle:", a3, "; of these also with (01):", both)
open("out_v5_s3.txt", "w").write("\n".join(out) + "\n")
