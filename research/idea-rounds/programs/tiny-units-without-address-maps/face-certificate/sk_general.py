# skeptic: ingredients of general.md recomputed with my own code (true matrices of regmachine.py).
import sys, itertools
import numpy as np
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..'))
import regmachine as R
np.seterr(all="ignore")
def wt(x): return bin(x).count("1")
def frank(M):
    if min(M.shape) == 0: return 0, float("inf")
    s = np.linalg.svd(M, compute_uv=False)
    if s[0] < 1e-9: return 0, float("inf")
    r = int((s > 1e-8*s[0]).sum()); gap = s[r-1]/s[r] if r < len(s) and s[r] > 0 else float("inf")
    return r, gap
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
def exact_crank(Z):
    A = np.rint(Z.real).astype(np.int64); B = np.rint(Z.imag).astype(np.int64)
    assert np.abs(Z - A - 1j*B).max() < 1e-6
    big = np.block([[A, -B], [B, A]]).tolist(); r = bareiss_rank(big); assert r % 2 == 0
    return r // 2
def minors(g, S, T, chunk=150000):
    S = np.array(S); T = np.array(T); out = np.zeros((len(S), len(T)), complex)
    step = max(1, chunk // max(1, len(T)))
    for a in range(0, len(S), step):
        sub = g[S[a:a+step, None, :, None], T[None, :, None, :]]
        out[a:a+step] = np.linalg.det(sub)
    return out
part = sys.argv[1]
if part == "thm2":
    print("THEOREM 2 of general.md: N_u = span of e(x_1..x_f) with u.x_1 = ... = u.x_f; whole space V, Kronecker powers of the true matrices")
    for m, fs in ((2, (3, 4, 5)), (3, (3, 4)), (4, (3,))):
        n = 1 << m
        for f in fs:
            def kp(g):
                M = g
                for _ in range(f-1): M = np.kron(M, g)
                return M
            tup = np.array(list(itertools.product(range(n), repeat=f)))
            K = kp(R.kernel(m)); bad = 0; lines = []
            dirs = {w: None for w in range(1, n)}
            for u in range(1, n):
                par = np.array([[R.dot(u, int(x)) for x in t] for t in tup]) if n**f <= 4096 else None
                pb = np.array([R.dot(u, x) for x in range(n)])[tup]
                inn = np.where((pb == pb[:, :1]).all(axis=1))[0]; out = np.where(~(pb == pb[:, :1]).all(axis=1))[0]
                dimN = len(inn); assert dimN == 2*(n//2)**f
                dem, gap = frank(K[np.ix_(out, inn)])
                pred = dimN if wt(u) % 2 else dimN - dimN // 2**(f-1)
                ok = dem == pred; bad += (not ok)
                lines.append((u, dimN, dem, pred, gap))
            # shifts and phase letters keep N_u; directions with u.w = 0 keep it; with u.w = 1: rank
            capbad = 0; keptbad = 0; capfull = 0; captested = 0; zero_tested = 0
            for w in range(1, n):
                Dw = kp(R.dirm(m, w))
                for u in range(1, n):
                    pb = np.array([R.dot(u, x) for x in range(n)])[tup]; same = (pb == pb[:, :1]).all(axis=1)
                    inn = np.where(same)[0]; out = np.where(~same)[0]
                    blk = Dw[np.ix_(out, inn)]
                    if R.dot(u, w) == 0:
                        zero_tested += 1; capbad += (np.abs(blk).max() > 1e-12)
                    elif (m < 4) or (u in (1, 3, 7, 15) and w in (1, 15, 6, 8)):
                        r, gap = frank(blk); captested += 1; capfull += (r == len(inn))
            free = [R.shift(m, s) for s in range(n)] + [R.phase_mat(m, z, c) for z in range(n) for c in (0, 1)]
            for g in free:
                G = kp(g)
                for u in (range(1, n) if m < 4 else (1, 3, 7, 15)):
                    pb = np.array([R.dot(u, x) for x in range(n)])[tup]; same = (pb == pb[:, :1]).all(axis=1)
                    keptbad += (np.abs(G[np.ix_(np.where(~same)[0], np.where(same)[0])]).max() > 1e-12)
            print(f"m={m} f={f}: dem as predicted for {len(lines)-bad} of {len(lines)} functionals u; least singular-value gap {min(l[4] for l in lines):.3g}; examples (u, dim N_u, dem, predicted): {[l[:4] for l in lines[:3]]}")
            print(f"      directions with u.w=0: block (outside N_u) x N_u non-zero in {capbad} of {zero_tested} cases; with u.w=1: cap = dim N_u in {capfull} of {captested} tested; shifts and phase letters moving N_u outside itself: {keptbad}")
            sys.stdout.flush()
if part == "bases":
    print("28 BASES at block size 3, three columns, antisymmetric part: N_U = span of a(T), T inside a plane u.x = c, u in U")
    m, n, f = 3, 8, 3
    sets = list(itertools.combinations(range(n), f))
    K8 = 8*R.kernel(m); D2 = {w: 2*R.dirm(m, w) for w in range(1, n)}
    bases = set()
    for U in itertools.combinations(range(1, n), 3):
        if U[0] ^ U[1] != U[2]: bases.add(U)
    print("number of bases of functionals:", len(bases))
    res = {}
    for U in sorted(bases):
        inn = [T for T in sets if any(len({R.dot(u, x) for x in T}) == 1 for u in U)]; out = [T for T in sets if T not in set(inn)]
        caps = []
        for w in range(1, n):
            M = minors(D2[w], out, inn); e = exact_crank(M); fl, gap = frank(M); assert e == fl; caps.append(e)
        M = minors(K8, out, inn); dem = exact_crank(M); fl, gap = frank(M); assert dem == fl
        res[U] = (len(inn), caps, dem)
    print("dim N_U values:", sorted({v[0] for v in res.values()}), "; largest cap over all bases and directions:", max(max(v[1]) for v in res.values()), "; dem values:", sorted({v[2] for v in res.values()}))
    print("coordinate basis (1,2,4):", res[(1, 2, 4)])
    print("bases with cap = 8 for all seven directions:", sum(1 for v in res.values() if v[1] == [8]*7), "of", len(res))
    print("bases with dem = 24 and every cap <= 8:", sum(1 for v in res.values() if v[2] == 24 and max(v[1]) <= 8), "of", len(res))
if part == "facet":
    print("FACET SUBSPACE, antisymmetric part: N = span of a(T), T an f-set inside a facet x_j = c (in at least one facet)")
    for m, f in ((2, 2), (3, 3), (3, 4), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7), (4, 8)):
        n = 1 << m
        sets = list(itertools.combinations(range(n), f))
        isin = [any(len({(x >> j) & 1 for x in T}) == 1 for j in range(m)) for T in sets]
        inn = [T for T, b in zip(sets, isin) if b]; out = [T for T, b in zip(sets, isin) if not b]
        byw = {}
        gapmin = float("inf")
        for w in range(1, n):
            if m == 4 and f >= 5 and w not in (1, 8, 3, 12, 7, 14, 15): continue
            r, gap = frank(minors(2*R.dirm(m, w), out, inn)); byw.setdefault(wt(w), set()).add(r); gapmin = min(gapmin, gap)
        dem, gap = frank(minors((2**m)*R.kernel(m), out, inn)); gapmin = min(gapmin, gap)
        keptbad = 0
        rng = np.random.default_rng(5)
        free = [R.shift(m, s) for s in range(n)] + [R.phase_mat(m, z, c) for z in range(n) for c in (0, 1)]
        for _ in range(6 if f >= 5 else 12):
            g = np.eye(n, dtype=complex)
            for i in rng.integers(0, len(free), 5): g = g @ free[i]
            keptbad += (np.abs(minors(g, out, inn)).max() > 1e-9)
        c = 2*len(list(itertools.combinations(range(n//2), f))) if f <= n//2 else 0
        print(f"m={m} f={f}: dim A {len(sets)}, dim N {len(inn)}, c = 2*C(2^(m-1),f) = {c}; cap by weight {dict((k, sorted(v)) for k, v in sorted(byw.items()))}; dem {dem}; random free products moving N: {keptbad}; least gap {gapmin:.3g}")
        sys.stdout.flush()
