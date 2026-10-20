# skeptic: independent recomputation of the four-column certificate at block size 4.
# Method: true matrices of regmachine.py, applied to the vectors b(P) as tensors of shape 16^4
# (one tensordot per column; no minors, no reduced model, no primes); exact Gaussian integers;
# rank = exact whole-number (Bareiss) rank of a Gram matrix.
import sys, json, itertools, random, re
from fractions import Fraction
import numpy as np
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..'))
import regmachine as R
np.seterr(all="ignore")
m, n, f = 4, 16, 4
sets = list(itertools.combinations(range(n), f))
flat = np.array([((T[0]*n+T[1])*n+T[2])*n+T[3] for T in sets])
def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i+1, len(p)):
            if p[i] > p[j]: s = -s
    return s
PERMS = [(p, perm_sign(p)) for p in itertools.permutations(range(f))]
def a_tensor(T):
    T = tuple(sorted(T)); v = np.zeros((n,)*f, complex)
    for p, s in PERMS:
        v[tuple(T[i] for i in p)] += s
    return v
def b_tensor(P):
    P = sorted(P); seen = set(); v = np.zeros((n,)*f, complex)
    for t in range(n):
        c = tuple(sorted(t ^ x for x in P))
        if c in seen: continue
        seen.add(c); v += a_tensor(c)
    assert len(seen) == 4
    return v
def apply(g, X):            # g on every column of each tensor X[r]
    for ax in range(1, f+1):
        X = np.moveaxis(np.tensordot(g, X, axes=(1, ax)), 0, ax)
    return X
def coords(X):              # coordinates on sorted tuples (injective on antisymmetric tensors)
    Y = X.reshape(X.shape[0], -1)[:, flat]
    Z = np.rint(Y.real) + 1j*np.rint(Y.imag)
    assert np.abs(Y - Z).max() < 1e-6, np.abs(Y - Z).max()
    return Z
def antisym_ok(X):          # is each tensor antisymmetric under exchange of columns 1,2 and 1..4 cyclic
    return (np.abs(X + np.swapaxes(X, 1, 2)).max() < 1e-9) and (np.abs(X + np.swapaxes(X, 3, 4)).max() < 1e-9) and (np.abs(X + np.swapaxes(X, 2, 3)).max() < 1e-9)
def bareiss_rank(M):
    M = [list(map(int, r)) for r in M]; nr = len(M); nc = len(M[0]) if nr else 0
    rank = 0; prev = 1
    for c in range(nc):
        piv = next((r for r in range(rank, nr) if M[r][c] != 0), None)
        if piv is None: continue
        M[rank], M[piv] = M[piv], M[rank]; p = M[rank][c]; Mp = M[rank]
        for r in range(rank+1, nr):
            q = M[r][c]; Mr = M[r]
            M[r] = [(p*Mr[j] - q*Mp[j]) // prev for j in range(nc)]
        prev = p; rank += 1
        if rank == nr: break
    return rank
def crank_exact(Z):         # exact rank over C of a Gaussian-integer matrix Z (rows = vectors)
    if Z.shape[0] == 0: return 0
    re_, im_ = np.rint(Z.real).astype(np.int64), np.rint(Z.imag).astype(np.int64)
    A = re_.astype(object); B = im_.astype(object)
    Gr = A.dot(A.T) + B.dot(B.T); Gi = B.dot(A.T) - A.dot(B.T)     # Gram = Z Z^H
    r = Z.shape[0]
    big = [[0]*(2*r) for _ in range(2*r)]
    for i in range(r):
        for j in range(r):
            big[i][j] = int(Gr[i, j]); big[i][j+r] = -int(Gi[i, j]); big[i+r][j] = int(Gi[i, j]); big[i+r][j+r] = int(Gr[i, j])
    rr = bareiss_rank(big); assert rr % 2 == 0
    return rr // 2
def crank_float(Z):
    if Z.shape[0] == 0: return 0, None
    s = np.linalg.svd(Z, compute_uv=False); r = int((s > 1e-9*s[0]).sum())
    gap = (s[r-1]/s[r] if r < len(s) and s[r] > 0 else float("inf")) if r > 0 else None
    return r, gap
def off_N(gB, B):           # 4*(gB) minus its part along N; the rows of B are orthogonal with norm^2 4
    nb = (np.abs(B)**2).sum(axis=1); assert np.all(nb == 4)
    return 4*gB - (gB @ B.conj().T) @ B
D2 = {w: 2*R.dirm(m, w) for w in range(1, n)}        # Gaussian integers
K16 = 16*R.kernel(m)
assert np.abs(K16 - np.rint(K16.real) - 1j*np.rint(K16.imag)).max() < 1e-12
P = np.eye(n)
for j in range(m): P = P @ R.dirm(m, 1 << j)
assert np.abs(P - R.kernel(m)).max() < 1e-12, "kernel is not the product of the unit directions"
FREE = [("shift", s, None, R.shift(m, s)) for s in range(n)] + [("phase", z, c, R.phase_mat(m, z, c)) for z in range(n) for c in (0, 1)]
rng = random.Random(20261011)
def random_free_product(L=6):
    g = np.eye(n, dtype=complex); names = []
    for _ in range(L):
        k, a, c, M = rng.choice(FREE); g = g @ M; names.append((k, a, c))
    return names, g
rows = []
txt = open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "cert_4col.txt")).read()
blocks = re.findall(r"row \[(.*?)\]: (\d+) planes: (\[\[.*?\]\])\n\s+cap for directions 1\.\.15: (\[.*?\]) ; summed by number of 1-bits: (\[.*?\]) ; dem (\d+) ; multiplier (\d+)/(\d+)", txt)
assert len(blocks) == 4
allok = True
out = []
for name, npl, planes, caps_a, sums_a, dem_a, mu_n, mu_d in blocks:
    planes = json.loads(planes); caps_a = json.loads(caps_a); sums_a = json.loads(sums_a); dem_a = int(dem_a)
    for Pl in planes:      # each is a 2-dimensional subspace
        assert Pl[0] == 0 and len(set(Pl)) == 4 and (Pl[1] ^ Pl[2]) == Pl[3]
    Bt = np.stack([b_tensor(Pl) for Pl in planes]); assert antisym_ok(Bt)
    B = coords(Bt); r = len(planes)
    assert crank_exact(B) == r
    # (1) kept: all 48 generators, and 40 random products of 6 generators
    bad = 0; tested = 0
    for k, a, c, M in FREE:
        gB = coords(apply(M, Bt)); tested += 1
        if np.abs(off_N(gB, B)).max() != 0: bad += 1
    for _ in range(40):
        names, g = random_free_product(); gB = coords(apply(g, Bt)); tested += 1
        if np.abs(off_N(gB, B)).max() != 0: bad += 1
    # control: an address map that is not a shift should move N for some
    G = R.perm_mat(m, R.linear_map(m, [1, 3, 4, 8]))
    ctl = np.abs(off_N(coords(apply(G, Bt)), B)).max() != 0
    # (2) ranks
    caps = []; gaps = []
    for w in range(1, n):
        gt = apply(D2[w], Bt); assert antisym_ok(gt)
        Rm = off_N(coords(gt), B); e = crank_exact(Rm); fl, gap = crank_float(Rm); assert e == fl, (w, e, fl)
        caps.append(e); gaps.append(gap)
    Rm = off_N(coords(apply(K16, Bt)), B); dem = crank_exact(Rm); fl, gapk = crank_float(Rm); assert dem == fl
    sums = [sum(caps[w-1] for w in range(1, n) if bin(w).count("1") == j) for j in (1, 2, 3, 4)]
    ok = (caps == caps_a and dem == dem_a and sums == sums_a)
    allok &= ok and bad == 0
    print(f"row [{name}]: planes {r}; kept: {bad} of {tested} free letters (48 generators + 40 random products of 6) move N outside itself; control address map moves N: {ctl}")
    print(f"   my cap 1..15 {caps}\n   author's     {caps_a}\n   my dem {dem} (author {dem_a}); my sums by weight {sums} (author {sums_a}); all equal: {ok}; least singular-value gap {min(g for g in gaps+[gapk] if g):.3g}")
    out.append({"name": name, "planes": planes, "caps": caps, "dem": dem, "sums": sums, "mult": [int(mu_n), int(mu_d)]})
json.dump(out, open("my_cert_rows.json", "w"))
# (3) certificate by plain sums
lam = [Fraction(o["mult"][0], o["mult"][1]) for o in out]
cnt = [4, 6, 4, 1]
tot = [sum(l*o["sums"][j] for l, o in zip(lam, out)) for j in range(4)]
rhs = sum(l*o["dem"] for l, o in zip(lam, out))
print("multipliers", [str(l) for l in lam], "all >= 0:", all(l >= 0 for l in lam))
print("weighted sum of my rows by weight:", [str(t) for t in tot], "; number of directions:", cnt, "; equal:", tot == cnt)
print("weighted sum of my demands:", rhs, "=", float(rhs), "; equals 3728/1219:", rhs == Fraction(3728, 1219), "; below 4:", rhs < 4)
print("ALL OK" if allok and tot == cnt and rhs == Fraction(3728, 1219) else "MISMATCH")
