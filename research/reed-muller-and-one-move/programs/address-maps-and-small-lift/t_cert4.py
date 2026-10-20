"""t_cert4.py (rm-theory): the block-4 certificate of wht38/tiny/close (wedge4/cert_4col.txt: four kept subspaces of
Lambda^4(C^16), multipliers 438/6095, 392/6095, 45/1219, 22/1219, bound 3728/1219 for the STANDARD kernel), evaluated on
the kernel in EVERY basis b of F_2^4 (420 unordered bases), as Theorem C3's transfer needs (NOTE 8.9, Next 3).

For each row r and basis b: dem_r(b) = dim(N_r + K_b N_r) - dim N_r, K_b = dir(b1) dir(b2) dir(b3) dir(b4), computed
modulo two primes p = 1 mod 4 (i -> sqrt(-1)); a rank modulo p is at most the true rank, so these are LOWER bounds on
the true demands, which is the safe direction. The caps are wedge4's exact per-direction values (from the file).
The certificate is symmetrised over the 24 coordinate permutations sigma (each sigma N_r is kept by shifts and
diagonal letters, cap_{sigma N}(w) = cap_N(sigma^-1 w), dem_{sigma N}(K_b) = dem_N(K_{sigma^-1 b})).
Output: D(b) = sum_r mult_r * avg_sigma dem_r(sigma^-1 b), the combined cap per direction (must be <= 1), and min_b D(b):
the no-phase game at block 4 then needs at least min_b D(b) moves per array for every address group.
Also re-checks wedge4's caps (mod p, lower bounds) and dem 16 for the standard basis as calibration.
"""
import itertools, sys, time
from fractions import Fraction
import numpy as np

n = 16
ROWS = [
 ([[0,1,4,5],[0,1,10,11],[0,1,12,13],[0,1,14,15],[0,2,5,7],[0,2,9,11],[0,2,12,14],[0,2,13,15],[0,3,12,15],[0,4,9,13],
   [0,4,10,14],[0,4,11,15],[0,5,9,12],[0,5,10,15],[0,5,11,14],[0,7,9,14],[0,7,10,13],[0,7,11,12]],
  [6,6,5,6,3,5,6,5,6,6,3,6,7,2,6], Fraction(438, 6095)),
 ([[0,1,2,3],[0,1,4,5],[0,1,6,7],[0,1,8,9],[0,2,4,6],[0,2,5,7],[0,2,9,11],[0,2,13,15],[0,3,4,7],[0,3,5,6],[0,3,8,11],
   [0,4,9,13],[0,4,11,15],[0,5,8,13],[0,6,9,15],[0,6,11,13],[0,7,8,15],[0,7,10,13]],
  [6,2,6,3,6,3,6,7,6,5,6,5,6,5,6], Fraction(392, 6095)),
 ([[0,1,6,7],[0,1,8,9],[0,1,12,13],[0,1,14,15],[0,2,12,14],[0,3,4,7],[0,3,8,11],[0,3,12,15],[0,3,13,14],[0,4,8,12],
   [0,4,9,13],[0,4,11,15],[0,6,8,14],[0,6,9,15],[0,6,10,12],[0,6,11,13],[0,7,8,15],[0,7,9,14],[0,7,11,12]],
  [6,5,6,6,5,6,3,3,6,5,6,6,7,6,3], Fraction(45, 1219)),
 ([[0,1,6,7],[0,1,14,15],[0,2,4,6],[0,2,13,15],[0,3,4,7],[0,3,5,6],[0,3,8,11],[0,3,9,10],[0,3,12,15],[0,4,8,12],
   [0,4,11,15],[0,5,10,15],[0,6,8,14],[0,6,9,15],[0,6,10,12],[0,6,11,13],[0,7,8,15]],
  [7,6,5,6,5,6,5,6,3,7,5,3,7,6,4], Fraction(22, 1219)),
]
QUADS = list(itertools.combinations(range(n), 4))
QIDX = {q: k for k, q in enumerate(QUADS)}
PERM4 = list(itertools.permutations(range(4)))
SGN = [(-1) ** sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j]) for p in PERM4]
QA = np.array(QUADS)                                   # 1820 x 4


def sqrtm1(p):
    for a in range(2, 10000):
        r = pow(a, (p - 1) // 4, p)
        if r * r % p == p - 1: return r


def rank_mod(A, p):
    A = A.copy() % p; r = 0; rows, cols = A.shape
    for c in range(cols):
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        piv = r + nz[0]
        if piv != r: A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): A[nzr] = (A[nzr] - np.outer(col[nzr], A[r])) % p
        r += 1
        if r == rows: break
    return r


def wedge(Mcols, p):
    """Mcols: 16 x 4 matrix (columns = images of the four basis vectors), mod p -> 1820 Plucker coordinates."""
    out = np.zeros(len(QUADS), dtype=np.int64)
    for pm, sg in zip(PERM4, SGN):
        t = np.ones(len(QUADS), dtype=np.int64)
        for k in range(4):
            t = (t * Mcols[QA[:, k], pm[k]]) % p
        out = (out + sg * t) % p
    return out


def plane_vec(P):
    """b(P) = sum over the four cosets t + P of a(t + P), with wedge4's convention (w4.py): a(t + P) is the shift by t
    of a(P) = e_p1 ^ e_p2 ^ e_p3 ^ e_p4 (p1 < ... < p4), i.e. eps * e_sorted with eps the sign of sorting."""
    Ps = sorted(P); v = {}; cos = []
    for t in range(n):
        img = [x ^ t for x in Ps]
        srt = tuple(sorted(img))
        inv = sum(1 for a in range(4) for b in range(a + 1, 4) if img[a] > img[b])
        sg = (-1) ** inv
        k = QIDX[srt]
        if k in v:
            assert v[k] == sg
            continue
        v[k] = sg; cos.append((srt, sg))
    return v, cos


def main():
    t0 = time.time()
    primes = [32749, 32717]
    bases = []
    for b in itertools.combinations(range(1, n), 4):
        basis = []
        for x in b:
            for bb in basis: x = min(x, x ^ bb)
            if x: basis.append(x)
        if len(basis) == 4: bases.append(b)
    print(f"{len(bases)} unordered bases of F_2^4", flush=True)
    results = {}
    for p in primes:
        I = sqrtm1(p)
        def dirm(w):
            M = np.zeros((n, n), dtype=np.int64)
            for x in range(n):
                M[x, x] = (M[x, x] + 1 + I) % p; M[x ^ w, x] = (M[x ^ w, x] + 1 - I) % p
            return M
        D = {w: dirm(w) for w in range(1, n)}
        rowdata = []
        for planes, caps, mult in ROWS:
            vecs = []; cosets = []
            for P in planes:
                v, cs = plane_vec(P); vecs.append(v); cosets.append(cs)
            Nmat = np.zeros((len(vecs), len(QUADS)), dtype=np.int64)
            for i, v in enumerate(vecs):
                for k, c in v.items(): Nmat[i, k] = c % p
            dimN = rank_mod(Nmat, p)
            rowdata.append((Nmat, cosets, dimN))
        def image(G, cosets):
            rows = []
            for cs in cosets:
                acc = np.zeros(len(QUADS), dtype=np.int64)
                for c, sg in cs:
                    acc = (acc + sg * wedge(G[:, list(c)], p)) % p
                rows.append(acc)
            return np.array(rows)
        def f_of(G, r):
            Nmat, cosets, dimN = rowdata[r]
            return rank_mod(np.concatenate([Nmat, image(G, cosets)]), p) - dimN
        # calibration: caps (mod p lower bounds) and the standard kernel
        capchk = [[f_of(D[w], r) for w in range(1, n)] for r in range(4)]
        Kstd = D[1] @ D[2] % p @ D[4] % p @ D[8] % p
        demstd = [f_of(Kstd, r) for r in range(4)]
        print(f"p={p}: dim N = {[rd[2] for rd in rowdata]}; caps mod p = {capchk}; match file: "
              f"{[capchk[r] == ROWS[r][1] for r in range(4)]}; dem(standard) = {demstd} ({time.time()-t0:.0f} s)",
              flush=True)
        dem = {}
        for bi, b in enumerate(bases):
            K = D[b[0]] @ D[b[1]] % p @ D[b[2]] % p @ D[b[3]] % p
            dem[b] = tuple(f_of(K, r) for r in range(4))
            if bi % 60 == 0: print(f"   basis {bi}/{len(bases)} {b}: {dem[b]} ({time.time()-t0:.0f} s)", flush=True)
        results[p] = dem
    dem = {b: tuple(min(results[p][b][r] for p in primes) for r in range(4)) for b in bases}
    # symmetrise over coordinate permutations
    def permute(w, s):
        y = 0
        for k in range(4):
            if (w >> k) & 1: y |= 1 << s[k]
        return y
    perms = list(itertools.permutations(range(4)))
    bkey = {tuple(sorted(b)): b for b in bases}
    Dval = {}
    for b in bases:
        tot = Fraction(0)
        for r in range(4):
            acc = 0
            for s in perms:
                bb = tuple(sorted(permute(w, s) for w in b))
                acc += dem[bkey[bb]][r]
            tot += ROWS[r][2] * Fraction(acc, len(perms))
        Dval[b] = tot
    capc = []
    for w in range(1, n):
        tot = Fraction(0)
        for r in range(4):
            tot += ROWS[r][2] * Fraction(sum(ROWS[r][1][permute(w, s) - 1] for s in perms), len(perms))
        capc.append(tot)
    vals = sorted(set(Dval.values()))
    print("combined cap per direction (must be <= 1):", sorted(set(capc)))
    print("combined demand over the 420 bases: distinct values", [str(v) for v in vals][:10], "min", vals[0],
          "=", float(vals[0]), "; standard basis:", Dval[(1, 2, 4, 8)], "=", float(Dval[(1, 2, 4, 8)]))
    hist = {}
    for b, v in Dval.items(): hist[str(v)] = hist.get(str(v), 0) + 1
    print("histogram:", hist)
    # GL-averaged certificate (Averaging Lemma over all address maps): every g N_r is kept by shifts and diagonals;
    # the direct sum over g in GL(4,2) has cap = average over the 15 directions, dem = average over the 840 bases.
    import json
    json.dump({",".join(map(str, b)): dem[b] for b in bases}, open("dem_cert4.json", "w"))
    avgdem = [Fraction(sum(dem[b][r] for b in bases), len(bases)) for r in range(4)]
    avgcap = [Fraction(sum(ROWS[r][1]), 15) for r in range(4)]
    comb = sum(ROWS[r][2] * avgdem[r] for r in range(4)) / sum(ROWS[r][2] * avgcap[r] for r in range(4))
    print("per row: average dem over the 840 bases", [str(x) for x in avgdem], "; average cap", [str(x) for x in avgcap])
    print("per row ratio avg dem / avg cap:", [float(avgdem[r] / avgcap[r]) for r in range(4)])
    print("GL-averaged combined certificate (wedge4 multipliers):", comb, "=", float(comb))
    print("best single GL-averaged row:", max(float(avgdem[r] / avgcap[r]) for r in range(4)))
    pairs_ok = all(dem[b1] == dem[b2] for b1 in bases for b2 in bases if False)
    print(f"done in {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()
