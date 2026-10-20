"""v_geo.py (v-w3-b5, skeptic): independent dual-polar-graph code for Sp(10,2), written from scratch (no sp2.py import).
Checks, over EVERY geodesic walk M0-L0-N1-N2 from one fixed adjacent pair (L0, M0) (no transitivity assumed):
  - anchors P: nbrs of M0 with d(P, N2) = 4 (and d(P, L0) = 2, d(P, N1) = 3 implied?)
  - ordered pairs (P1, Q): A1 in nbr(N2) with d(P1,A1) = 5, d(Q,A1) <= 4; A2 symmetric; census of lambda layers
  - Lemma L3's exact fact: for every pair admitting lambda = (3,3), every Y in nbr(N2) with d(M0, Y) = 4 has
    max(d(P1,Y), d(Q,Y)) = 5.
Also tabulates the same quantity for pairs admitting lambda_1 = 3 only (for robustness), and the general list.
Sp(2m,2) is transitive on ordered adjacent pairs, so fixing (L0, M0) loses nothing.
"""
import sys, time
from collections import Counter

m = 5
MX = (1 << m) - 1


def bc(x):
    return bin(x).count('1')


def swap(v):
    return ((v & MX) << m) | (v >> m)


def om(u, v):
    return bin(u & swap(v)).count('1') & 1


def span(basis):
    els = [0]
    for b in basis:
        els = els + [e ^ b for e in els]
    return els


def mask(basis):
    s = 0
    for e in span(basis):
        s |= 1 << e
    return s


def solve_dual(basis):
    """d_j with om(b_i, d_j) = delta_ij.  om(b, d) = parity(swap(b) & d)... careful: om(u,v)=par(u & swap(v)) = par(swap(u) & v)."""
    rows = [swap(b) for b in basis]
    n = len(rows)
    # augmented: rows with rhs identity
    aug = [[rows[i], 1 << i] for i in range(n)]
    piv_cols = []
    r = 0
    for c in range(2 * m):
        pr = next((k for k in range(r, n) if (aug[k][0] >> c) & 1), None)
        if pr is None:
            continue
        aug[r], aug[pr] = aug[pr], aug[r]
        for k in range(n):
            if k != r and (aug[k][0] >> c) & 1:
                aug[k][0] ^= aug[r][0]; aug[k][1] ^= aug[r][1]
        piv_cols.append(c); r += 1
        if r == n:
            break
    assert r == n
    duals = []
    for j in range(n):
        d = 0
        for k in range(n):
            if (aug[k][1] >> j) & 1:
                d |= 1 << piv_cols[k]
        duals.append(d)
    for i in range(n):
        for j in range(n):
            assert om(basis[i], duals[j]) == (1 if i == j else 0)
    return duals


class Lag:
    __slots__ = ('basis', 'mask')

    def __init__(self, basis):
        self.basis = tuple(basis)
        self.mask = mask(basis)


def dist(A, B):
    k = bc(A.mask & B.mask)
    return m - (k.bit_length() - 1)


_nbr_cache = {}


def nbrs(X):
    key = X.mask
    if key in _nbr_cache:
        return _nbr_cache[key]
    b = X.basis
    d = solve_dual(list(b))
    out = []
    for f in range(1, 1 << m):
        # basis of U_f = ker f on X: vectors sum c_i b_i with c.f = 0
        i0 = (f & -f).bit_length() - 1
        U = []
        for i in range(m):
            if i == i0:
                continue
            if (f >> i) & 1:
                U.append(b[i] ^ b[i0])
            else:
                U.append(b[i])
        df = 0
        for j in range(m):
            if (f >> j) & 1:
                df ^= d[j]
        for w in (df, df ^ b[i0]):
            Y = Lag(U + [w])
            out.append(Y)
    _nbr_cache[key] = out
    return out


def check_lag(X):
    for u in X.basis:
        for v in X.basis:
            assert om(u, v) == 0
    assert (bc(X.mask)) == 32


def main():
    t0 = time.time()
    L0 = Lag([1 << i for i in range(m)])        # pure-X Lagrangian (any Lagrangian will do)
    check_lag(L0)
    nL0 = nbrs(L0)
    assert len({Y.mask for Y in nL0}) == 62 and all(dist(L0, Y) == 1 for Y in nL0)
    M0 = nL0[int(sys.argv[1]) if len(sys.argv) > 1 else 7]
    nM0 = nbrs(M0)
    walks = []
    for N1 in nL0:
        if dist(M0, N1) != 2:
            continue
        for N2 in nbrs(N1):
            if dist(L0, N2) == 2 and dist(M0, N2) == 3:
                walks.append((N1, N2))
    N2set = {}
    for N1, N2 in walks:
        N2set.setdefault(N2.mask, N2)
    print("walks", len(walks), "distinct N2", len(N2set), flush=True)
    census = Counter(); factc = Counter(); fact1 = Counter(); anchc = Counter(); pairc = Counter(); lamc = Counter()
    npairs33 = 0; viol = 0; viol1 = 0
    farc = Counter(); impl = Counter()
    for idx, N2 in enumerate(N2set.values()):
        nN2 = nbrs(N2)
        assert len({Y.mask for Y in nN2}) == 62
        anchors = [P for P in nM0 if dist(P, N2) == 4]
        anchc[len(anchors)] += 1
        for P in anchors:
            impl[(dist(P, L0),)] += 1
        dM = [dist(M0, Y) for Y in nN2]
        far = 0
        for j, x in enumerate(dM):
            if x == 4:
                far |= 1 << j
        farc[bc(far)] += 1
        F5 = []; F3 = []; F4 = []; Lf = []
        for P in anchors:
            f5 = f3 = f4 = 0
            for j, Y in enumerate(nN2):
                x = dist(P, Y)
                if x == 5:
                    f5 |= 1 << j
                elif x == 3:
                    f3 |= 1 << j
                elif x == 4:
                    f4 |= 1 << j
            F5.append(f5); F3.append(f3); F4.append(f4)
            Lf.append((f3 | f4) & far)       # far Y with d(P, Y) <= 4 (d >= 3 always since d(P,N2)=4)
        n = len(anchors)
        np_ = 0
        for i in range(n):
            for k in range(n):
                if i == k:
                    continue
                a1_3 = F5[i] & F3[k]; a1_4 = F5[i] & F4[k]
                a2_3 = F5[k] & F3[i]; a2_4 = F5[k] & F4[i]
                c1 = a1_3 | a1_4; c2 = a2_3 | a2_4
                if c1 and c2:
                    np_ += 1
                    lamc[(3, 3)] += bc(a1_3) * bc(a2_3)
                    lamc[(3, 4)] += bc(a1_3) * bc(a2_4)
                    lamc[(4, 3)] += bc(a1_4) * bc(a2_3)
                    lamc[(4, 4)] += bc(a1_4) * bc(a2_4)
                    pairc[(bc(c1), bc(c2))] += 1
                if a1_3 and a2_3:
                    npairs33 += 1
                    bad = Lf[i] & Lf[k]
                    if bad:
                        viol += 1
                if a1_3:                       # lambda_1 = 3 only (array 2 may be lambda 4)
                    if Lf[i] & Lf[k]:
                        viol1 += 1
                    else:
                        fact1['ok'] += 1
        census[np_] += 1
        if idx % 100 == 0:
            print(idx, "N2 done, t=%.1fs" % (time.time() - t0), flush=True)
    print("anchors per N2:", dict(anchc), " d(P,L0) of anchors:", dict(impl))
    print("far Y count per N2:", dict(farc))
    print("special ordered pairs per N2:", dict(census))
    print("per special pair (#clean A1, #clean A2):", dict(pairc))
    print("lambda configurations total (over all N2, x3 for walks):", dict(lamc))
    print("pairs admitting (3,3):", npairs33, " violations of the L3 fact:", viol)
    print("pairs admitting lambda_1=3 (any lambda_2): violations", viol1, "ok", fact1['ok'])
    print("time %.1fs" % (time.time() - t0))


if __name__ == '__main__':
    main()
