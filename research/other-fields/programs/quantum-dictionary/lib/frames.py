#!/usr/bin/env python3
"""frames.py - the finite frame graphs of the tiny frame game, their symmetries and pattern classes.

Conventions.  A vector of F_2^(2m) is an int v = x | (z << m) and stands for the Pauli X^x Z^z (up to phase);
omega(u, v) = x.z' + x'.z.  The label of a frame matrix S is L(S) = { P : S P S^-1 is diagonal } = S^-1 L_0 S with
L_0 = {(0|z)} (the diagonal Paulis).  dir(z0) acts on labels as the transvection along (z0|0), the kernel as
F : (x|z) -> (x+z|z), free matrices fix L_0.

  G_lag(m): all Lagrangian subspaces of F_2^(2m); edge iff dim(L & L') = m-1; session iff L & L' = 0.
  G_sub(m): all subspaces U of F_2^m;  edge iff nested with codimension 1; session iff U (+) V = F_2^m.
            U sits in G_lag(m) as L_U = {(u | u+v) : u in U, v in U^perp}.
Run "python3 frames.py" for the checks (also part of tests.py).
"""
import os
import sys
import itertools
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
if os.path.isdir(os.path.join(HERE, "pylib")) and os.path.join(HERE, "pylib") not in sys.path:
    sys.path.insert(0, os.path.join(HERE, "pylib"))


def popc(x):
    return bin(x).count("1")


def span(vs):
    s = {0}
    for v in vs:
        if v not in s:
            s |= {x ^ v for x in s}
    return frozenset(s)


def basis_of(U):
    b, sp = [], {0}
    for v in sorted(U):
        if v not in sp:
            b.append(v)
            sp |= {x ^ v for x in sp}
    return b


def dim2(U):
    return len(U).bit_length() - 1


def omega(u, v, m):
    M = (1 << m) - 1
    return (popc((u & M) & (v >> m)) + popc((v & M) & (u >> m))) & 1


def linmap(N, f):
    img = [f(1 << b) for b in range(N)]
    T = [0] * (1 << N)
    for v in range(1, 1 << N):
        low = v & -v
        T[v] = T[v ^ low] ^ img[low.bit_length() - 1]
    return T


def sp_generators(m):
    """(name, table on F_2^(2m), fixes L_0).  gl: address maps x_i += x_j; ph: phases i^Q; F: the kernel."""
    M = (1 << m) - 1
    bit = lambda v, i: (v >> i) & 1
    gens = []
    for i in range(m):
        for j in range(m):
            if i != j:
                def f(v, i=i, j=j):
                    x, z = v & M, v >> m
                    return (x ^ (bit(x, j) << i)) | ((z ^ (bit(z, i) << j)) << m)
                gens.append(("gl%d%d" % (i, j), linmap(2 * m, f), True))
    for i in range(m):
        for j in range(i, m):
            def f(v, i=i, j=j):
                x, z = v & M, v >> m
                z ^= (bit(x, i) << i) if i == j else ((bit(x, j) << i) ^ (bit(x, i) << j))
                return x | (z << m)
            gens.append(("ph%d%d" % (i, j), linmap(2 * m, f), True))
    gens.append(("F", linmap(2 * m, lambda v: ((v & M) ^ (v >> m)) | ((v >> m) << m)), False))
    return gens


_cache = {}


def lag_graph(m):
    if ("lag", m) in _cache:
        return _cache[("lag", m)]
    gens = sp_generators(m)
    L0 = frozenset(z << m for z in range(1 << m))
    LF = frozenset(z | (z << m) for z in range(1 << m))
    seen, fr = {L0}, [L0]
    while fr:
        nx = []
        for L in fr:
            for _, T, _ in gens:
                Lg = frozenset(T[v] for v in L)
                if Lg not in seen:
                    seen.add(Lg)
                    nx.append(Lg)
        fr = nx
    frames = sorted(seen, key=lambda s: sorted(s))
    idx = {L: i for i, L in enumerate(frames)}
    n = len(frames)
    masks = [sum(1 << v for v in L) for L in frames]
    adj, opp, half = [[] for _ in range(n)], [[] for _ in range(n)], 1 << (m - 1)
    for i in range(n):
        mi = masks[i]
        for j in range(n):
            c = popc(mi & masks[j])
            if c == half:
                adj[i].append(j)
            elif c == 1:
                opp[i].append(j)
    perms = [(nm, [idx[frozenset(T[v] for v in L)] for L in frames], fx) for nm, T, fx in gens]
    G = dict(kind="lag", m=m, n=n, frames=frames, idx=idx, adj=adj, opp=opp, perms=perms, i0=idx[L0], iF=idx[LF],
             desc=[basis_of(L) for L in frames], transitive=True)
    _cache[("lag", m)] = G
    return G


def perp(U, m):
    return frozenset(v for v in range(1 << m) if all(popc(u & v) % 2 == 0 for u in U))


def L_of(U, m):
    return frozenset(u | ((u ^ v) << m) for u in U for v in perp(U, m))


def sub_graph(m):
    if ("sub", m) in _cache:
        return _cache[("sub", m)]
    found, fr = {span([])}, [span([])]
    while fr:
        nx = []
        for U in fr:
            for v in range(1, 1 << m):
                if v not in U:
                    V = span(list(U) + [v])
                    if V not in found:
                        found.add(V)
                        nx.append(V)
        fr = nx
    frames = sorted(found, key=lambda s: (len(s), sorted(s)))
    idx = {U: i for i, U in enumerate(frames)}
    n = len(frames)
    adj = [[j for j in range(n) if (frames[i] < frames[j] or frames[j] < frames[i])
            and max(len(frames[i]), len(frames[j])) == 2 * min(len(frames[i]), len(frames[j]))] for i in range(n)]
    opp = [[j for j in range(n) if len(frames[i] & frames[j]) == 1 and len(frames[i]) * len(frames[j]) == 1 << m]
           for i in range(n)]
    perms = []
    for i in range(m):
        for j in range(m):
            if i != j:
                T = linmap(m, lambda v, i=i, j=j: v ^ (((v >> j) & 1) << i))
                perms.append(("gl%d%d" % (i, j), [idx[frozenset(T[v] for v in U)] for U in frames], False))
    perms.append(("perp", [idx[perp(U, m)] for U in frames], False))
    G = dict(kind="sub", m=m, n=n, frames=frames, idx=idx, adj=adj, opp=opp, perms=perms, desc=[basis_of(U) for U in frames],
             transitive=False)
    _cache[("sub", m)] = G
    return G


def graph(kind, m):
    return lag_graph(m) if kind == "lag" else sub_graph(m)


def bfs_dist(adj, s):
    d = {s: 0}
    fr = [s]
    while fr:
        nx = []
        for u in fr:
            for v in adj[u]:
                if v not in d:
                    d[v] = d[u] + 1
                    nx.append(v)
        fr = nx
    return d


# ----------------------------------------------------------------------------------------------------------------
# pattern classes: orbits of multisets of W start frames under the group generated by G["perms"]
# ----------------------------------------------------------------------------------------------------------------
def _components(n, k, perms):
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    N = comb(n + k - 1, k)
    arr = np.fromiter(itertools.chain.from_iterable(itertools.combinations_with_replacement(range(n), k)),
                      dtype=np.int64, count=N * k).reshape(N, k)
    pw = (n ** np.arange(k - 1, -1, -1)).astype(np.int64)
    code = arr @ pw
    rows, cols = [], []
    for p in perms:
        img = np.sort(np.asarray(p, dtype=np.int64)[arr], axis=1) @ pw
        rows.append(np.arange(N))
        cols.append(np.searchsorted(code, img))
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    nc, lab = connected_components(coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(N, N)), directed=False)
    return arr, code, pw, nc, lab


def pattern_classes(G, W, method="auto"):
    """sorted list of (representative tuple, orbit size or None)"""
    import numpy as np
    n = G["n"]
    perms = [p for _, p, _ in G["perms"]]
    if method == "auto":
        method = "direct" if comb(n + W - 1, W) <= 3000000 or not G["transitive"] else "pinned"
    if method == "direct":
        arr, code, pw, nc, lab = _components(n, W, perms)
        _, first, counts = np.unique(lab, return_index=True, return_counts=True)
        return sorted((tuple(int(x) for x in arr[i]), int(c)) for i, c in zip(first, counts))
    # pinned: every orbit meets the multisets containing i0; orbits of the rest under the stabiliser of i0, then
    # merged through elements carrying another member of the multiset to i0
    i0 = G["i0"]
    stab = [p for _, p, fx in G["perms"] if fx]
    if W == 1:
        return [((i0,), n)]
    arr, code, pw, nc, lab = _components(n, W - 1, stab)
    tr = {i0: list(range(n))}                      # tr[x][x] = i0
    fr = [i0]
    while fr:
        nx = []
        for x in fr:
            for p in perms:
                y = p[x]
                if y not in tr:
                    pinv = [0] * n
                    for a, b in enumerate(p):
                        pinv[b] = a
                    tr[y] = [tr[x][pinv[v]] for v in range(n)]
                    nx.append(y)
        fr = nx
    assert len(tr) == n
    parent = list(range(nc))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    _, first = np.unique(lab, return_index=True)
    for i in first:
        ms = [int(x) for x in arr[i]]
        for x in set(ms):
            if x == i0:
                continue
            img = sorted(tr[x][v] for v in ms + [i0])
            img.remove(i0)
            j = int(np.searchsorted(code, int(np.dot(img, pw))))
            a, b = find(int(lab[i])), find(int(lab[j]))
            if a != b:
                parent[max(a, b)] = min(a, b)
    reps = {}
    for i in first:
        r = find(int(lab[i]))
        t = tuple(sorted([int(x) for x in arr[i]] + [i0]))
        if r not in reps or t < reps[r]:
            reps[r] = t
    return sorted((t, None) for t in reps.values())


# ----------------------------------------------------------------------------------------------------------------
# m = 2: explicit 4x4 matrices over the Gaussian rationals.  A matrix is stored as a Gaussian-integer matrix
# (entries (re, im), Python ints) that is a power of two times the true one; every test below is scale-free.
# ----------------------------------------------------------------------------------------------------------------
def m_mul(A, B):
    k = len(A)
    out = []
    for i in range(k):
        row = []
        for j in range(k):
            re = im = 0
            for t in range(k):
                a, b = A[i][t], B[t][j]
                if a[0] or a[1]:
                    re += a[0] * b[0] - a[1] * b[1]
                    im += a[0] * b[1] + a[1] * b[0]
            row.append((re, im))
        out.append(row)
    return out


def m_dag(A):
    return [[(A[j][i][0], -A[j][i][1]) for j in range(len(A))] for i in range(len(A))]


Z0, O1 = (0, 0), (1, 0)
M_id = [[O1 if i == j else Z0 for j in range(4)] for i in range(4)]
M_shift = lambda z: [[O1 if (i ^ z) == j else Z0 for j in range(4)] for i in range(4)]            # BinaryFrames.shift
M_dir2 = lambda z: [[(1, 1) if i == j else ((1, -1) if (i ^ z) == j else Z0) for j in range(4)] for i in range(4)]  # 2*dir(z)
M_perm = lambda pi: [[O1 if pi[x] == y else Z0 for y in range(4)] for x in range(4)]              # Letters.permMat
M_phase = lambda z, c: [[((0, 1) if (popc(z & x) + c) & 1 else O1) if x == y else Z0 for y in range(4)] for x in range(4)]
M_Z = lambda a: [[((-1, 0) if popc(a & x) & 1 else O1) if x == y else Z0 for y in range(4)] for x in range(4)]


def M_kernel4():                                                                                  # 4 * littleC (x) littleC
    def ent(x, y):
        r = (1, 0)
        for b in range(2):
            c = (1, 1) if ((x >> b) & 1) == ((y >> b) & 1) else (1, -1)
            r = (r[0] * c[0] - r[1] * c[1], r[0] * c[1] + r[1] * c[0])
        return r
    return [[ent(x, y) for y in range(4)] for x in range(4)]


def pauli_of(Q):
    """Q = c * X^x Z^z  ->  x | z << 2, else None"""
    nz = [j for j in range(4) if Q[0][j] != Z0]
    if len(nz) != 1:
        return None
    x, c = nz[0], Q[0][nz[0]]
    sg = []
    for i in range(4):
        for j in range(4):
            if j != (i ^ x):
                if Q[i][j] != Z0:
                    return None
            elif Q[i][j] == c:
                sg.append(0)
            elif Q[i][j] == (-c[0], -c[1]):
                sg.append(1)
            else:
                return None
    if sg[3] != sg[1] ^ sg[2]:
        return None
    return x | ((sg[1] | (sg[2] << 1)) << 2)


def label_of(S):
    """the Lagrangian { P : S P S^-1 diagonal } = span of S^dagger Z^a S (S is unitary up to scale)"""
    Sd = m_dag(S)
    vs = [pauli_of(m_mul(m_mul(Sd, M_Z(a)), S)) for a in (1, 2)]
    assert None not in vs, "not a Clifford matrix"
    L = span(vs)
    assert len(L) == 4 and omega(vs[0], vs[1], 2) == 0
    return L


def to_compact(M):
    """monomial matrix with entries t * i^e (one positive t) -> (permutation, exponents), else None"""
    pi, ex, t = [], [], None
    for row in M:
        nz = [j for j in range(4) if row[j] != Z0]
        if len(nz) != 1:
            return None
        re, im = row[nz[0]]
        if re and im:
            return None
        a = abs(re) + abs(im)
        if t is None:
            t = a
        if a != t:
            return None
        pi.append(nz[0])
        ex.append(0 if re > 0 else 2 if re < 0 else 1 if im > 0 else 3)
    return tuple(pi), tuple(ex)


def free_group_m2():
    """closure of the generators of Free (Work/GFrame/Engine/Free.lean) at m = 2, as compact monomial matrices"""
    gens = [M_perm([0, 2, 1, 3]), M_perm([0, 1, 3, 2]), M_phase(1, 0), M_phase(2, 0), M_phase(3, 0), M_phase(0, 1),
            M_shift(1), M_shift(2)]
    cg = [to_compact(g) for g in gens]
    mul = lambda a, b: (tuple(b[0][a[0][x]] for x in range(4)), tuple((a[1][x] + b[1][a[0][x]]) % 4 for x in range(4)))
    e = (tuple(range(4)), (0, 0, 0, 0))
    seen, fr = {e}, [e]
    while fr:
        nx = []
        for a in fr:
            for g in cg:
                b = mul(a, g)
                if b not in seen:
                    seen.add(b)
                    nx.append(b)
        fr = nx
    return seen, gens


def check_m2_matrices():
    """Every edge of G_lag(2) is one dir(z) factor up to free matrices; equal labels differ by a free matrix;
    the kernel sends a frame to a transversal one and all 8 transversal frames occur.  Returns counts."""
    G = lag_graph(2)
    free, fgens = free_group_m2()
    assert len(free) == 3072, len(free)
    for g in fgens:                                    # free matrices keep the label L_0, hence every label
        assert label_of(g) == G["frames"][G["i0"]]
    ph = [M_id]
    for z in (1, 2, 3):
        ph += [m_mul(M_phase(z, 0), A) for A in ph]    # the 8 phases i^Q modulo Paulis
    lin = [[0, 1, 2, 3], [0, 2, 1, 3], [0, 1, 3, 2], [0, 3, 1, 2], [0, 2, 3, 1], [0, 3, 2, 1]]      # GL(2,2) on addresses
    reps = {G["i0"]: M_id}
    order, edges, coinc, prods = [G["i0"]], set(), 0, 0
    for a in order:
        S = reps[a]
        for z in (1, 2, 3):
            for A in ph:
                N = m_mul(M_dir2(z), m_mul(A, S))
                prods += 1
                b = G["idx"][label_of(N)]
                assert b in G["adj"][a], "a dir factor did not move to an adjacent frame"
                edges.add((a, b))
                if b not in reps:
                    reps[b] = N
                    order.append(b)
                else:                                   # equal labels: the ratio must be a free matrix (up to scale)
                    R = to_compact(m_mul(N, m_dag(reps[b])))
                    assert R is not None and R in free, "equal labels but the ratio is not free"
                    coinc += 1
    assert len(reps) == 15 and len(edges) == 90 == sum(len(x) for x in G["adj"])
    for a in range(15):                                 # free generators keep every label
        for g in fgens:
            assert G["idx"][label_of(m_mul(g, reps[a]))] == a
    K4, sess = M_kernel4(), 0
    for a in range(15):
        tg = set()
        for A in ph:
            for pi in lin:
                tg.add(G["idx"][label_of(m_mul(K4, m_mul(M_perm(pi), m_mul(A, reps[a]))))])
        assert tg == set(G["opp"][a]), "kernel targets are not exactly the transversal frames"
        sess += len(tg)
    # dir(z)^2 = shift(z):  (2 dir)^2 = 4 shift
    for z in (1, 2, 3):
        assert m_mul(M_dir2(z), M_dir2(z)) == [[(4 * e[0], 4 * e[1]) for e in r] for r in M_shift(z)]
    return dict(free_group=len(free), products=prods, directed_edges=len(edges), coincidences=coinc, sessions=sess)


# ----------------------------------------------------------------------------------------------------------------
def check_graphs(ms_lag=(2, 3, 4), ms_sub=(2, 3, 4), out=print):
    res = {}
    for m in ms_lag:
        G = lag_graph(m)
        n, adj, opp = G["n"], G["adj"], G["opp"]
        want = 1
        for i in range(1, m + 1):
            want *= 2 ** i + 1
        assert n == want, (n, want)
        for nm, T, _ in sp_generators(m):               # generators are symplectic
            assert all(omega(T[1 << a], T[1 << b], m) == omega(1 << a, 1 << b, m) for a in range(2 * m) for b in range(2 * m))
        for nm, p, fx in G["perms"]:                    # and act as automorphisms preserving edges and sessions
            assert sorted(p) == list(range(n)) and (p[G["i0"]] == G["i0"]) == fx
            for i in range(n):
                assert sorted(p[j] for j in adj[i]) == adj[p[i]] and sorted(p[j] for j in opp[i]) == opp[p[i]]
        assert all(len(a) == 2 * (2 ** m - 1) for a in adj) and all(len(o) == 2 ** (m * (m + 1) // 2) for o in opp)
        # one dir factor = the transvection along (z|0); the kernel = the product of m of them = F
        i0, L0 = G["i0"], G["frames"][G["i0"]]
        tau = lambda v, z: v ^ (omega(v, z, m) * z)
        assert G["idx"][frozenset(tau(v, 1) for v in L0)] in adj[i0]
        Ftab = [T for nm, T, _ in sp_generators(m) if nm == "F"][0]
        for v in range(1 << (2 * m)):
            w = v
            for i in range(m):
                w = tau(w, 1 << i)
            assert w == Ftab[v]
        # the images of the free matrices (the generators fixing L_0) are transitive on the neighbours of L_0 and on
        # the frames transversal to L_0: every edge is free * dir * free, every transversal pair is a session
        stab = [p for _, p, fx in G["perms"] if fx]

        def orbit(x):
            seen, fr = {x}, [x]
            while fr:
                y = fr.pop()
                for p in stab:
                    if p[y] not in seen:
                        seen.add(p[y])
                        fr.append(p[y])
            return seen
        assert orbit(adj[i0][0]) == set(adj[i0]) and orbit(G["iF"]) == set(opp[i0])
        srcs = range(n) if m <= 3 else [G["i0"]]
        for s in srcs:                                  # graph distance = m - dim(L & L'); sessions at distance m
            d = bfs_dist(adj, s)
            assert len(d) == n
            for t in range(n):
                assert d[t] == m - dim2(G["frames"][s] & G["frames"][t])
            assert all(d[t] == m for t in opp[s]) and [t for t in range(n) if d[t] == m] == opp[s]
        npairs = None
        if m <= 3:                                      # sessions (g L_0, g L_F) are exactly the transversal pairs
            seen, fr = {(G["i0"], G["iF"])}, [(G["i0"], G["iF"])]
            while fr:
                nx = []
                for a, b in fr:
                    for _, p, _ in G["perms"]:
                        q = (p[a], p[b])
                        if q not in seen:
                            seen.add(q)
                            nx.append(q)
                fr = nx
            assert seen == {(a, b) for a in range(n) for b in opp[a]}
            npairs = len(seen)
        order = None
        if m == 2:                                      # the group used for symmetry reduction is Sp(4,2), order 720
            ps = [tuple(p) for _, p, _ in G["perms"]]
            seen, fr = {tuple(range(n))}, [tuple(range(n))]
            while fr:
                nx = []
                for a in fr:
                    for g in ps:
                        b = tuple(g[x] for x in a)
                        if b not in seen:
                            seen.add(b)
                            nx.append(b)
                fr = nx
            order = len(seen)
            assert order == 720
        res[("lag", m)] = dict(frames=n, degree=len(adj[0]), sessions_per_frame=len(opp[0]), session_pairs=npairs,
                               generators=len(G["perms"]), group_order=order, distance_sources=len(srcs))
        out("  G_lag(%d): %s" % (m, res[("lag", m)]))
    for m in ms_sub:
        G = sub_graph(m)
        n, adj, opp, fr_ = G["n"], G["adj"], G["opp"], G["frames"]
        assert n == {2: 5, 3: 16, 4: 67}[m]
        for nm, p, _ in G["perms"]:
            assert sorted(p) == list(range(n))
            for i in range(n):
                assert sorted(p[j] for j in adj[i]) == adj[p[i]] and sorted(p[j] for j in opp[i]) == opp[p[i]]
        Lg = lag_graph(m)
        emb = [Lg["idx"][L_of(U, m)] for U in fr_]      # L_U is one of the Lagrangian frames
        assert len(set(emb)) == n
        F = [p for nm, p, _ in Lg["perms"] if nm == "F"][0]
        ns = 0
        for i in range(n):
            assert F[emb[i]] == emb[G["idx"][perp(fr_[i], m)]]          # kernel: L_U -> L_(U perp)
            d = bfs_dist(adj, i)
            for j in range(n):
                dl = m - dim2(Lg["frames"][emb[i]] & Lg["frames"][emb[j]])
                assert d[j] == dl == dim2(fr_[i]) + dim2(fr_[j]) - 2 * dim2(fr_[i] & fr_[j])
                assert (j in adj[i]) == (emb[j] in Lg["adj"][emb[i]]) and (j in opp[i]) == (emb[j] in Lg["opp"][emb[i]])
                assert (j in opp[i]) == (d[j] == m)
            ns += len(opp[i])
        res[("sub", m)] = dict(frames=n, edges=sum(len(a) for a in adj) // 2, session_pairs=ns, generators=len(G["perms"]))
        out("  G_sub(%d): %s" % (m, res[("sub", m)]))
    return res


if __name__ == "__main__":
    import time
    t = time.time()
    print("m = 2 explicit matrices:", check_m2_matrices())
    check_graphs()
    for kind, m, W in (("sub", 2, 3), ("lag", 2, 2), ("lag", 2, 3), ("lag", 2, 4), ("sub", 3, 3), ("lag", 3, 2), ("lag", 3, 3)):
        print("  classes %s m=%d W=%d: %d" % (kind, m, W, len(pattern_classes(graph(kind, m), W))))
    a = pattern_classes(lag_graph(2), 4, "pinned"), pattern_classes(lag_graph(3), 3, "pinned")
    print("  pinned method: lag m=2 W=4: %d, lag m=3 W=3: %d" % (len(a[0]), len(a[1])))
    print("all graph checks passed in %.1fs" % (time.time() - t))
