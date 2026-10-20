# vc_ab.py -- v-cert-07: independent rebuild of cert-07's intersection-size cuts N_I = span{e_S (x) e_T : |S cap T| in I}
# in X = Lam^a V (x) Lam^b V (block 4, V = C^16), for a + b ODD, by a GLOBAL reduction (not cert-07's per-block sectors):
# the Pauli centre acts by (-1)^(a+b) = -1, so X = V (x) M' with M' = span{[S,T] : sig S = sig T} (Z-character 0),
# Phi(e_y (x) [S,T]) = (Lam^a X^y (x) Lam^b X^y)(e_S (x) e_T).  For G = Lam^a g (x) Lam^b g (g Clifford):
# rho(G)[(S',T'),(S,T)] = det g[S',S] det g[T',T] / (lam g00)   (only y = 0 can hit a Z-character-0 basis tensor since
# sig(S+y) + sig(T+y) = (a+b) y = y).  N_I = V (x) (N_I cap M'), so psi_j(full) = 16 * rank in M'.
# Depth element: real rotation r = [[1,1],[-1,1]]/sqrt2 on address bits 0..j-1 (diagonal conjugate of the engine's d).
# usage: vc_ab.py a b I(comma list)  [j list, default 1,3]
import sys, itertools, time, json
import numpy as np, os
t0 = time.time()
MB = int(os.environ.get('VC_M', '4')); NV = 1 << MB   # block size m (default 4)
a, b = int(sys.argv[1]), int(sys.argv[2]); I = set(int(x) for x in sys.argv[3].split(','))
JS = [int(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else [1, 3]
assert (a + b) % 2 == 1
def sig(S):
    r = 0
    for s in S: r ^= s
    return r
SA_ = list(itertools.combinations(range(NV), a)); TB_ = list(itertools.combinations(range(NV), b))
byS = {}
for T in TB_: byS.setdefault(sig(T), []).append(T)
PAIRS = [(S, T) for S in SA_ for T in byS.get(sig(S), [])]
NP = len(PAIRS)
inN = np.array([len(set(S) & set(T)) in I for S, T in PAIRS])
print(f"(a,b,I) = ({a},{b},{sorted(I)}): dim M' = {NP} (expect C(n,a)C(n,b)/n = {len(SA_)*len(TB_)//NV}), dim N cap M' = {int(inN.sum())}", flush=True)
r = np.array([[1.0, 1.0], [-1.0, 1.0]]) / np.sqrt(2)
res = {}
for j in JS:
    n = 1 << j; Z = 1 << (MB - j)
    D = np.array([[1.0]])
    for _ in range(j): D = np.kron(D, r)
    pc = [bin(c).count('1') for c in range(1 << n)]
    Mn = np.zeros((1 << n, 1 << n))
    for c1 in range(1 << n):
        rr = [i for i in range(n) if c1 >> i & 1]
        for c2 in range(1 << n):
            if pc[c2] != pc[c1]: continue
            cc = [i for i in range(n) if c2 >> i & 1]
            Mn[c1, c2] = np.linalg.det(D[np.ix_(rr, cc)]) if rr else 1.0
    def codes(sets):
        out = np.zeros((len(sets), Z), dtype=np.int64)
        for k, S in enumerate(sets):
            for s in S: out[k, s >> j] |= 1 << (s & (n - 1))
        return out
    cS = codes([p[0] for p in PAIRS]); cT = codes([p[1] for p in PAIRS])
    key = [tuple(pc[x] for x in cS[k]) + tuple(pc[x] for x in cT[k]) for k in range(NP)]
    cls = {}
    for k, kk in enumerate(key): cls.setdefault(kk, []).append(k)
    psi = 0; phi = 0.0; lo = 9.0; hi = 0.0; uerr = 0.0; big = 0
    for kk, rows in cls.items():
        rows = np.array(rows); m = inN[rows]
        big = max(big, len(rows))
        if m.all() or (~m).all():
            continue
        src = rows[m]; dst = rows[~m]
        G = np.ones((len(dst), len(src)))
        for z in range(Z):
            G *= Mn[cS[dst, z][:, None], cS[src, z][None, :]]
            G *= Mn[cT[dst, z][:, None], cT[src, z][None, :]]
        G /= D[0, 0]
        s = np.linalg.svd(G, compute_uv=False)
        tol = 1e-7
        psi += int((s > tol).sum()); phi += float((s ** 2).sum())
        if (s > tol).any(): lo = min(lo, float(s[s > tol].min()))
        if (s <= tol).any(): hi = max(hi, float(s[s <= tol].max()))
        if len(rows) <= 2000:   # unitarity spot check of the whole class block
            Gf = np.ones((len(rows), len(rows)))
            for z in range(Z):
                Gf *= Mn[cS[rows, z][:, None], cS[rows, z][None, :]]
                Gf *= Mn[cT[rows, z][:, None], cT[rows, z][None, :]]
            Gf /= D[0, 0]
            uerr = max(uerr, float(np.abs(Gf @ Gf.T - np.eye(len(rows))).max()))
    res[j] = (psi, phi)
    print(f"  j={j}: {len(cls)} classes (max {big}), psi = {psi} (x{NV} = {NV*psi}), phi = {phi:.6f} (x{NV} = {NV*phi:.3f}), "
          f"s = {phi/psi if psi else 0:.6f}, kept >= {lo:.3e}, dropped <= {hi:.1e}, unitarity err (classes <= 2000) {uerr:.1e} ({time.time()-t0:.0f}s)", flush=True)
if 1 in res and 3 in res:
    print(f"RESULT ({a},{b},{sorted(I)}): m = {MB}: psi1 = {NV*res[1][0]}, psi3 = {NV*res[3][0]}, ratio3 = {res[3][0]/res[1][0]:.6f}, "
          f"spec3 = {res[3][1]/res[1][1]:.6f}", flush=True)
