# hookrel.py -- cert-07 (a): hook relatives in Lam^a V (x) Lam^b V, V = C^16 (block 4), N = span{e_S (x) e_T : for every
# element e of [16], the pattern (e in S, e in T) lies in ALLOWED}  (hook: a = 1, allowed = all but (1,0), i.e. S subset T).
# Such N is invariant under every monomial map of V (all free steps of the engine and widened games).
# Reduction [proof]: g_j = h^{(x)j} (x) I_{2^(4-j)} splits V = (+)_z V_z (Z = 2^(4-j) blocks of size n' = 2^j), g_j acts
# by H_j = h^{(x)j} on every block, and Lam^a V (x) Lam^b V = (+)_{(a_z, b_z)} (x)_z (Lam^{a_z} V_z (x) Lam^{b_z} V_z), each
# summand preserved by G = Lam^a g (x) Lam^b g.  Inside a summand N is the tensor product of the per-block pieces F_z
# (the condition is per element).  For tensor products, dim(N + GN) - dim N = prod f_z - prod (f_z - rho_z) and
# ||(1-P)GP||_F^2 = prod f_z - prod (f_z - phi_z) (dim(A1(x)A2 cap B1(x)B2) = prod dim(Ai cap Bi); ||PGP||^2 multiplies).
# So psi_j(a, b) = [x^a y^b] ( F(x,y)^Z - (F - R)(x,y)^Z ), F = sum f(a',b') x^a' y^b', R = sum rho(a',b') x^a' y^b',
# rho, phi = rank and Frobenius mass of the block-j KERNEL H_j on the small piece (n' = 2, 4, 8: dense, exact-ish SVD).
import itertools, sys, json, time
import numpy as np
from fractions import Fraction as Fr
from math import comb
def lam(Hm, a):
    n = Hm.shape[0]; subs = list(itertools.combinations(range(n), a))
    if a == 0: return np.ones((1, 1)), subs
    M = np.empty((len(subs), len(subs)))
    for i, S in enumerate(subs):
        rows = Hm[list(S)]
        for k, T in enumerate(subs): M[i, k] = np.linalg.det(rows[:, list(T)])
    return M, subs
h = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
PATS = [(0, 0), (0, 1), (1, 0), (1, 1)]
def small_tables(j, allowed):
    n = 2 ** j; Hj = np.array([[1.0]])
    for _ in range(j): Hj = np.kron(Hj, h)
    L = [lam(Hj, a) for a in range(n + 1)]
    f = {}; rho = {}; phi = {}; gap = {}
    for a in range(n + 1):
        for b in range(n + 1):
            A, SA = L[a]; B, SB = L[b]
            inF = np.array([all(((e in S), (e in T)) in allowed for e in range(n)) for S in SA for T in SB])
            f[a, b] = int(inF.sum())
            if f[a, b] == 0 or f[a, b] == len(inF): rho[a, b] = 0; phi[a, b] = 0.0; continue
            G = np.kron(A, B)
            M = G[np.ix_(~inF, inF)]
            s = np.linalg.svd(M, compute_uv=False)
            r = int((s > 1e-9).sum()); rho[a, b] = r; phi[a, b] = float((s ** 2).sum())
            gap[a, b] = (float(s[r - 1]) if r else None, float(s[r]) if r < len(s) else None)
    return f, rho, phi, gap
def polypow_coef(T, Z, n, exact=True):
    """T: dict (a',b') -> coefficient; returns dict (a,b) -> coefficient of T(x,y)^Z, degrees <= 16"""
    P = {(0, 0): (1 if exact else 1.0)}
    for _ in range(Z):
        Q = {}
        for (a1, b1), c1 in P.items():
            for (a2, b2), c2 in T.items():
                if c2 == 0: continue
                k = (a1 + a2, b1 + b2)
                if k[0] > 16 or k[1] > 16: continue
                Q[k] = Q.get(k, 0) + c1 * c2
        P = Q
    return P
def profile_all(allowed, js=(1, 2, 3)):
    out = {}
    for j in js:
        f, rho, phi, gap = small_tables(j, allowed)
        Z = 2 ** (4 - j)
        PF = polypow_coef(f, Z, 2 ** j)
        PR = polypow_coef({k: f[k] - rho[k] for k in f}, Z, 2 ** j)
        PP = polypow_coef({k: f[k] - phi[k] for k in f}, Z, 2 ** j, exact=False)
        mg = min((g[0] for g in gap.values() if g[0]), default=None)
        out[j] = dict(psi={k: PF.get(k, 0) - PR.get(k, 0) for k in PF}, phi={k: PF.get(k, 0) - PP.get(k, 0.0) for k in PF},
                      dim={k: PF[k] for k in PF}, mingap=mg)
    return out
if __name__ == '__main__':
    t0 = time.time()
    # calibration: hook a = 1, b = 8 (27456, 48576, 68320) and b = 7 (24024, ?, 59078)
    allowed = set(PATS) - {(1, 0)}
    pr = profile_all(allowed)
    for b in (7, 8):
        print("hook k =", b, "psi_1..3 =", [pr[j]['psi'][(1, b)] for j in (1, 2, 3)], "phi_1..3 =",
              [round(pr[j]['phi'][(1, b)], 3) for j in (1, 2, 3)], "dimN", pr[1]['dim'][(1, b)], flush=True)
    rows = []
    for fb in range(16):
        allowed = {p for i, p in enumerate(PATS) if not (fb >> i) & 1}
        if len(allowed) in (0, 4): continue
        pr = profile_all(allowed)
        for (a, b), d in sorted(pr[1]['dim'].items()):
            if a > 16 or b > 16 or d == 0: continue
            tot = comb(16, a) * comb(16, b)
            if d == tot: continue
            p1, p3 = pr[1]['psi'].get((a, b), 0), pr[3]['psi'].get((a, b), 0)
            if p1 == 0: continue
            f1, f3 = pr[1]['phi'][(a, b)], pr[3]['phi'][(a, b)]
            rows.append(dict(forbid=[PATS[i] for i in range(4) if (fb >> i) & 1], a=a, b=b, dimN=d, dimX=tot,
                             psi=[p1, pr[2]['psi'][(a, b)], p3], phi=[f1, pr[2]['phi'][(a, b)], f3],
                             ratio3=p3 / p1, spec3=f3 / f1 if f1 else None, s1=f1 / p1, s3=f3 / p3 if p3 else None,
                             mingap=[pr[j]['mingap'] for j in (1, 2, 3)]))
        print("forbid", [PATS[i] for i in range(4) if (fb >> i) & 1], "done", f"{time.time()-t0:.0f}s", flush=True)
    rows.sort(key=lambda r: -r['ratio3'])
    with open("out_hookrel.txt", "w") as fh:
        for r in rows[:60]:
            line = (f"forbid {r['forbid']} a={r['a']} b={r['b']} dimN {r['dimN']}/{r['dimX']}: psi1..3 {r['psi']} ratio3 "
                    f"{Fr(r['psi'][2], r['psi'][0])} = {r['ratio3']:.5f} spec3 {r['spec3']:.5f} s1 {r['s1']:.5f} s3 {r['s3']:.5f}")
            fh.write(line + "\n"); print(line)
    json.dump(rows, open("out_hookrel.json", "w"), default=str)
    print(f"total {time.time()-t0:.0f}s, {len(rows)} certificates")
