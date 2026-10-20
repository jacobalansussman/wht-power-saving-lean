#!/usr/bin/env python3
"""cutconst.py (sl-combinatorics): the coherent-class cut constant c_m of the frame graph Gamma_m = dual polar graph
DSp(2m,2) (Lagrangians of F_2^(2m), adjacent iff they meet in dimension m-1).

PROOF.md Theorem 6: every coherent unit at block m (any number of scratch arrays) has M >= c_m W, where
    c_m = max over vertex sets Q of  [share of opposite pairs (distance m) split by Q] / [share of edges cut by Q].
If c_m = m, no coherent unit saves at block m.  c_m >= 2(1 - 2^-m) always (the sets {E : w in E}).
Exhaustive over all 2^15 sets at m = 2; structured families (balls, Voronoi half-spaces, vector sets) and a local
search at m = 3, 4.   usage: python3 cutconst.py m [local-search restarts]"""
import sys, random, itertools, time
import numpy as np

def lagrangians(m):
    n = 2 * m
    def omega(u, v):
        a, b = u & ((1 << m) - 1), u >> m
        c, d = v & ((1 << m) - 1), v >> m
        return (bin(a & d).count("1") + bin(b & c).count("1")) & 1
    def span_set(vecs):
        S = {0}
        for v in vecs: S |= {s ^ v for s in S}
        return S
    def bits(S):
        x = 0
        for s in S: x |= 1 << s
        return x
    L0 = span_set([1 << (m + j) for j in range(m)])            # Z labels (0|b)
    start = bits(L0)
    seen = {start}; frontier = [L0]
    vecs = list(range(1, 1 << n))
    while frontier:
        nf = []
        for L in frontier:
            for v in vecs:
                img = {u ^ (v if omega(u, v) else 0) for u in L}   # transvection T_v
                b = bits(img)
                if b not in seen: seen.add(b); nf.append(img)
        frontier = nf
    return sorted(seen), start

def setup(m):
    labs, L0 = lagrangians(m)
    V = len(labs)
    A = np.array([[bin(x & y).count("1") for y in labs] for x in labs], dtype=np.int64)
    dim = np.round(np.log2(A)).astype(np.int64)
    D = m - dim
    return labs, L0, D

def ratio(D, m, Q):
    q = np.zeros(D.shape[0], dtype=bool); q[list(Q)] = True
    split = q[:, None] != q[None, :]
    e = (D == 1); o = (D == m)
    fe = (split & e).sum() / e.sum(); fo = (split & o).sum() / o.sum()
    return (fo / fe if fe > 0 else 0.0), fo, fe

def main():
    m = int(sys.argv[1]); restarts = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    t0 = time.time()
    labs, L0, D = setup(m)
    V = len(labs)
    print("m=%d: |Gamma| = %d, valency %d, opposite per vertex %d, diameter %d (%.1f s)" %
          (m, V, (D[0] == 1).sum(), (D[0] == m).sum(), D.max(), time.time() - t0), flush=True)
    best = (0, None, "")
    def consider(Q, name):
        nonlocal best
        r, fo, fe = ratio(D, m, Q)
        if r > best[0] + 1e-12: best = (r, sorted(Q), name)
        return r
    # balls around a vertex
    for k in range(m):
        Q = np.flatnonzero(D[0] <= k); consider(Q, "ball radius %d" % k)
        print("  ball radius %d: ratio %.6f" % (k, ratio(D, m, Q)[0]))
    # Voronoi half-spaces between two vertices at distance j
    for j in range(1, m + 1):
        F2 = int(np.flatnonzero(D[0] == j)[0])
        for tie in ("<", "<="):
            Q = np.flatnonzero(D[0] < D[F2]) if tie == "<" else np.flatnonzero(D[0] <= D[F2])
            r = consider(Q, "Voronoi d(.,F1) %s d(.,F2), d(F1,F2)=%d" % (tie, j))
            print("  Voronoi %s, d(F1,F2)=%d: ratio %.6f" % (tie, j, r))
    # vector sets {E : w in E} (the colouring family) and intersections of two
    labs_arr = labs
    for w in (1, 1 << m, 1 | (1 << m)):
        Q = [i for i, x in enumerate(labs_arr) if (x >> w) & 1]
        r = consider(Q, "vector set w=%d" % w); print("  vector set w=%d: ratio %.6f (2(1-2^-m) = %.6f)" % (w, r, 2 * (1 - 2.0 ** -m)))
    if m == 2:
        # exhaustive over all subsets
        e = (D == 1); o = (D == m)
        best_ex = 0; arg = None
        for mask in range(1, 1 << (V - 1)):
            q = np.array([(mask >> i) & 1 for i in range(V)], dtype=bool)
            split = q[:, None] != q[None, :]
            fe = (split & e).sum(); fo = (split & o).sum()
            if fe == 0: continue
            r = (fo / o.sum()) / (fe / e.sum())
            if r > best_ex + 1e-12: best_ex = r; arg = np.flatnonzero(q).tolist()
        print("  EXHAUSTIVE over all %d subsets: c_2 = %.6f, attained by %s" % ((1 << (V - 1)) - 1, best_ex, arg))
        consider(arg, "exhaustive optimum")
    else:
        rng = random.Random(1)
        e = (D == 1); o = (D == m)
        ne, no = e.sum(), o.sum()
        for rs in range(restarts):
            q = np.array([rng.random() < 0.5 for _ in range(V)], dtype=bool)
            if rs < len(best[1] or []) or rs == 0:
                if best[1] is not None and rs == 0:
                    q = np.zeros(V, dtype=bool); q[best[1]] = True
            # cut counts per vertex flip, incremental
            improved = True
            def score(q):
                split = q[:, None] != q[None, :]
                fe = (split & e).sum(); fo = (split & o).sum()
                return ((fo / no) / (fe / ne)) if fe else 0.0
            cur = score(q); it = 0
            while improved and it < 60:
                improved = False; it += 1
                order = list(range(V)); rng.shuffle(order)
                # flip gains: for vertex v, edges to same side become cut and vice versa
                same = (q[:, None] == q[None, :])
                ce = ((~same) & e).sum() // 2; co = ((~same) & o).sum() // 2
                de = np.where(same, 1, -1) * e; do = np.where(same, 1, -1) * o
                gain_e = de.sum(axis=1); gain_o = do.sum(axis=1)
                ne_new = ce + gain_e; no_new = co + gain_o
                with np.errstate(divide="ignore", invalid="ignore"):
                    sc = np.where(ne_new > 0, (no_new / (no / 2)) / (ne_new / (ne / 2)), 0)
                v = int(np.argmax(sc))
                if sc[v] > cur + 1e-12:
                    q[v] = ~q[v]; cur = sc[v]; improved = True
            consider(np.flatnonzero(q), "local search restart %d" % rs)
        print("  local search (%d restarts): best %.6f" % (restarts, best[0]))
    print("RESULT m=%d: c_m >= %.6f (m = %d; plain needs %d; engine floor %s) from: %s" %
          (m, best[0], m, m, {2: "2", 3: "2.8", 4: "3.1538"}.get(m, "?"), best[2]))
    print("  best set size %d of %d" % (len(best[1]), V))

if __name__ == "__main__":
    main()
