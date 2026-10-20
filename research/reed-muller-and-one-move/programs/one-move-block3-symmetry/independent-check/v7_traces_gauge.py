import itertools, numpy as np
from vlib import *
exec(open("v2_f21.py").read().split("# all elements g")[0].split("c2 = c @ c")[1])  # dq, SYMS
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
L = lagrangians(); d1 = [j for j in range(len(L)) if len(L[j] & LZ) == 4]
def bitperm(p):
    return permm(lambda x: sum(((x >> j) & 1) << p[j] for j in range(3)))
t3 = bitperm((1, 2, 0)); t2 = bitperm((1, 0, 2))
def reps_for(cells):
    reps = {}
    for Sm in SYMS:
        D = dq(Sm)
        for A in gl3():
            PA = permm(lambda x, A=A: lin(A, x))
            for z in range(1, 8):
                m = D @ PA @ dirm(z); cl = apply_sym(sym(m), LZ)
                if cl in cells and cl not in reps: reps[cl] = m
            if len(reps) == len(cells): return reps
    return reps
for nm, g in [("3-cycle", t3), ("(01)", t2)]:
    cells = [L[j] for j in d1 if apply_sym(sym(g), L[j]) == L[j]]
    reps = reps_for(cells)
    for cl, m in reps.items():
        NM = orth(np.hstack([Nb, rho(m) @ Nb]))
        tr = quot_trace(rho(g), NM, Nb)
        tr3 = quot_trace(rho(t3), NM, Nb) if apply_sym(sym(t3), cl) == cl else None
        tr2 = quot_trace(rho(t2), NM, Nb) if apply_sym(sym(t2), cl) == cl else None
        say(nm, "-fixed cell", sorted(cl)[:4], "...: trace of 3-cycle", None if tr3 is None else np.round(tr3, 6), " trace of (01)", None if tr2 is None else np.round(tr2, 6))
# gauge-symmetric single move, f = 1 and f = 2: Z7 regular orbit of arrays, c Singer, c' = K^-1 c^-1 K
c = np.load("v_c.npy"); cp = Kinv @ inv(c) @ Kmat
rng = np.random.default_rng(15)
m = dirm(3)
for f in (1, 2):
    def tens(g): return g if f == 1 else np.kron(g, g)
    y = [rng.normal(size=8 ** f) + 1j * rng.normal(size=8 ** f) for _ in range(7)]
    # register value as coefficient list a_l (group-algebra elements = matrices here); x = sum_l c'^l e_l
    a = [np.linalg.matrix_power(cp, l) for l in range(7)]
    val = lambda coeffs, yy: sum(tens(coeffs[l]) @ yy[l] for l in range(7))
    mx = [m @ a[l] for l in range(7)]
    worst = 0
    for h in range(1, 7):
        ch = np.linalg.matrix_power(c, h); cph = np.linalg.matrix_power(cp, h)
        # tau_h maps coefficient a_l on array l to c_h a_l c'_h on array h+l
        tau = [None] * 7
        for l in range(7): tau[(l + h) % 7] = ch @ mx[l] @ cph
        worst = max(worst, np.abs(val(tau, y) - tens(ch) @ val(mx, y)).max())
    say("f = %d: max |tau_h(m x) - c_h (m x)| over h = 1..6 on random raw arrays: %.2e" % (f, worst))
open("out_v7_traces_gauge.txt", "w").write("\n".join(out) + "\n")
