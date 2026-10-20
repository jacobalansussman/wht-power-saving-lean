import itertools, numpy as np
from vlib import *
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
c = np.load("v_c.npy"); EN = np.load("v_EN.npy")
L = lagrangians(); d1 = [j for j in range(len(L)) if len(L[j] & LZ) == 4]
c2 = c @ c
def dq(Sm):
    d = []
    for x in range(8):
        v = [(x >> j) & 1 for j in range(3)]
        q = sum(Sm[j][j] * v[j] for j in range(3)) + 2 * sum(Sm[j][k] * v[j] * v[k] for j in range(3) for k in range(j + 1, 3))
        d.append(1j ** (q % 4))
    return np.diag(d)
SYMS = [[[b[0], b[3], b[4]], [b[3], b[1], b[5]], [b[4], b[5], b[2]]] for b in itertools.product([0, 1], repeat=6)]
def fixscalar(g, order):
    h = np.linalg.matrix_power(g, order); lam = h[0, 0]
    if not np.allclose(h, lam * np.eye(8)): return None
    for k in range(order * 8):
        t = np.exp(2j * np.pi * k / (order * 8)) * abs(lam) ** (-1 / order)
        if np.allclose(np.linalg.matrix_power(g * t, order), np.eye(8)): return g * t
    return None
# all elements g = P_A D_S X^a Z^b (B mod scalars) that are in the Levi and satisfy g c g^-1 = c^2 exactly (up to scalar 1)
cands = []
for A in gl3():
    PA = permm(lambda x, A=A: lin(A, x))
    for Sm in SYMS:
        base = PA @ dq(Sm)
        for a in range(8):
            for b in range(8):
                g = base @ pau(a, b)
                r = g @ c @ inv(g)
                if np.allclose(r, c2):
                    if is_free(Kinv @ g @ Kmat):
                        g3 = fixscalar(g, 3)
                        if g3 is not None: cands.append(g3)
say("Levi elements n (mod scalar) with n c n^-1 = c^2 exactly and n^3 = 1:", len(cands))
n = cands[0]
cp = Kinv @ inv(c) @ Kmat; cpn = Kinv @ inv(n) @ Kmat
assert is_free(cpn)
Rcpn = [rho(np.linalg.matrix_power(cpn, s)) for s in range(3)]
d = mults([quot_trace(R, Nb, EN) for R in Rcpn], 3); say("d = split of c'_n on N/E_N:", d)
Rc = rho(cp); P0 = sum(np.linalg.matrix_power(Rc, h) for h in range(7)) / 7
N0 = orth(P0 @ Nb); E0 = orth(P0 @ EN)
p = mults([quot_trace(R, N0, E0) for R in Rcpn], 3); say("p = split of c'_n on c'-invariant part of N/E_N (dims", N0.shape[1], E0.shape[1], "):", p)
# n-fixed cells at distance 1
sn = sym(n)
fix1 = [j for j in d1 if apply_sym(sn, L[j]) == L[j]]
say("n fixes", len(fix1), "of the 14 distance-1 cells")
# find letters in each fixed cell: m = g dir1 with g free (permutation * phase)
def cell_of(m): return apply_sym(sym(m), LZ)
reps = {}
for A in gl3():
    PA = permm(lambda x, A=A: lin(A, x))
    for z in range(1, 8):
        m = PA @ dirm(z)
        cl = cell_of(m)
        if cl in [L[j] for j in fix1] and cl not in reps: reps[cl] = m
    if len(reps) == len(fix1): break
say("representative letters found for", len(reps), "fixed cells")
Rn = [rho(np.linalg.matrix_power(n, s)) for s in range(3)]
def sigma_of(gm, m):
    NM = orth(np.hstack([Nb, rho(m) @ Nb]))
    Rg = [rho(np.linalg.matrix_power(gm, s)) for s in range(3)]
    # check stability
    assert np.linalg.matrix_rank(np.hstack([NM, Rg[1] @ NM]), tol=1e-7) == NM.shape[1]
    return mults([quot_trace(R, NM, Nb) for R in Rg], 3)
for cl, m in reps.items():
    say("fixed cell", sorted(cl), ": sigma(n) on (N+mN)/N =", sigma_of(n, m))
# exact commutation: search beta in B with m beta commuting with n
DS = np.array([dq(Sm) @ pau(a, b) for Sm in SYMS for a in range(8) for b in range(8)])
tot = {sorted(cl)[1]: [0, 0] for cl in reps}
nB = 0
for A in gl3():
    PA = permm(lambda x, A=A: lin(A, x))
    chunk = PA[None] @ DS; nB += len(chunk)
    for cl, m in reps.items():
        mm = m[None] @ chunk
        lhs = n[None] @ mm; rhs = mm @ n[None]
        fl = lhs.reshape(len(mm), -1); fr = rhs.reshape(len(mm), -1)
        k0 = np.argmax(np.abs(fr), axis=1)
        ratio = fl[np.arange(len(mm)), k0] / fr[np.arange(len(mm)), k0]
        ok = np.all(np.abs(fl - ratio[:, None] * fr) < 1e-8, axis=1)
        ex = ok & (np.abs(ratio - 1) < 1e-8)
        tot[sorted(cl)[1]][0] += int(ok.sum()); tot[sorted(cl)[1]][1] += int(ex.sum())
say("|B mod scalars| scanned:", nB)
for k, v in tot.items(): say("cell with second element", k, ": letters in the coset commuting with n up to scalar:", v[0], " exactly:", v[1])
# Pauli-twisted order-3 lifts of the same linear map (gauge alternatives): g = n P with (nP)^3 = 1 (up to scalar)
tw = []
for a in range(8):
    for b in range(8):
        g = fixscalar(n @ pau(a, b), 3)
        if g is not None: tw.append(((a, b), g))
say("Pauli twists n X^a Z^b of order 3 (mod scalar):", len(tw))
splits = {}
for (ab, g) in tw:
    dd = mults([quot_trace(rho(np.linalg.matrix_power(Kinv @ inv(g) @ Kmat, s)), Nb, EN) for s in range(3)], 3)
    ss = tuple(tuple(sigma_of(g, m)) for cl, m in reps.items())
    splits.setdefault((tuple(dd), ss), []).append(ab)
for k, v in splits.items(): say("twist demand split", k[0], " cap splits on the fixed cells", k[1], " count", len(v))
np.save("v_n.npy", n)
open("out_v2_f21.txt", "w").write("\n".join(out) + "\n")
