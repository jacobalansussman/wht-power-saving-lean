#!/usr/bin/env python3
"""trsub.py (w8-rmexact): Theorem R over GL(m,2)-submodules of the permutation blocks of C[Sym_m x| GL(m,2)].
Block O (a GL-orbit of binary quadratic forms): V_O = C[O], move [z z^T] = diag((-1)^{Q(z)}), target [I] =
diag((-1)^{tr Q}), GL permutes O.  Y_0 = an H-submodule (H = GL), f(g) = rank((1-P) Pi_g P) (Pi_g = the -1 projector).
Theorem R: M f(move) >= W f(target).  We decompose C[O] into irreducible H-submodules (eigenspaces of a generic
Hermitian element of the commutant), enumerate all sums of them (exact when C[O] is multiplicity free), and report
the best ratio f(I)/f(zz^T).  Output: per block, constituents (dim), multiplicity-free flag, best ratio."""
import sys, itertools, numpy as np
sys.path.insert(0, '.')
from forms import setup
m = int(sys.argv[1]) if len(sys.argv) > 1 else 4
D = setup(m); GL = D['GL']; supp = D['supp']; nf = D['nf']; ap = D['apply']
key = {supp[Q]: Q for Q in range(nf)}
dirs = D['dirs']
# permutation of forms by each h: Q -> form with support h(supp Q)
def act(h, Q): return key[frozenset(ap(h, x) for x in supp[Q])]
Itr = [sum(D['Qval'](Q, 1 << i) for i in range(m)) & 1 for Q in range(nf)]
rng = np.random.default_rng(1)
seen = [False] * nf
best_all = 0
for Q0 in range(1, nf):
    if seen[Q0]: continue
    # orbit and coset reps
    rep = {Q0: None}; order = [Q0]
    perms = []
    for h in GL:
        Q2 = act(h, Q0)
        if Q2 not in rep: rep[Q2] = h; order.append(Q2)
    for Q2 in order: seen[Q2] = True
    n = len(order); pos = {Q: i for i, Q in enumerate(order)}
    stab = [h for h in GL if act(h, Q0) == Q0]
    # suborbits of stab on the orbit
    sub = [-1] * n; ns = 0
    for i, Q in enumerate(order):
        if sub[i] >= 0: continue
        for h in stab: sub[pos[act(h, Q)]] = ns
        ns += 1
    # inverse of a GL element
    def inv(h):
        for g in GL:
            if all(ap(g, ap(h, 1 << i)) == 1 << i for i in range(m)): return g
    # label matrix L[a][b] = suborbit of t_a^-1 b
    L = np.zeros((n, n), dtype=int)
    for a, Qa in enumerate(order):
        t = rep[Qa]; ti = inv(t) if t is not None else None
        for b, Qb in enumerate(order):
            L[a, b] = sub[pos[act(ti, Qb)]] if ti is not None else sub[b]
    r = rng.normal(size=ns) + 1j * rng.normal(size=ns)
    T = r[L]; T = T + T.conj().T
    w, V = np.linalg.eigh(T)
    # cluster eigenvalues
    groups = []; cur = [0]
    for i in range(1, n):
        if abs(w[i] - w[i - 1]) < 1e-7 * max(1, abs(w[i])): cur.append(i)
        else: groups.append(cur); cur = [i]
    groups.append(cur)
    dims = [len(g) for g in groups]
    # multiplicity-freeness: commutant dim (= ns) equals number of constituents iff multiplicity free
    mf = (ns == len(groups))
    Piz = np.diag([float(D['Qval'](Q, 1)) for Q in order])
    PiI = np.diag([float(Itr[Q]) for Q in order])
    Bfull = np.hstack([V[:, g] for g in groups])
    Cz = Bfull.conj().T @ Piz @ Bfull; CI = Bfull.conj().T @ PiI @ Bfull
    cols = []; o = 0
    for g in groups: cols.append(list(range(o, o + len(g)))); o += len(g)
    def rk(X):
        if X.size == 0: return 0
        s = np.linalg.svd(X, compute_uv=False)
        return int((s > 1e-8).sum())
    best = (0, None, None, None)
    k = len(groups)
    if k > 16:
        print("block |O|=%d: %d constituents, too many to enumerate" % (n, k)); continue
    for mask in range(1, (1 << k) - 1):
        J = [c for j in range(k) if (mask >> j) & 1 for c in cols[j]]
        Jc = [c for j in range(k) if not (mask >> j) & 1 for c in cols[j]]
        fz = rk(Cz[np.ix_(Jc, J)]); fI = rk(CI[np.ix_(Jc, J)])
        if fz == 0 and fI > 0:
            best = (float('inf'), mask, fI, fz); continue
        if fz and fI / fz > best[0]: best = (fI / fz, mask, fI, fz)
    # isotypic components: group eigenspaces by dimension (for GL(4,2) = A_8 the repeated degrees 14, 20, 56 are unique)
    from collections import defaultdict
    iso = defaultdict(list)
    for j, g in enumerate(groups): iso[len(g)].append(j)
    Tp = (rng.normal(size=ns) + 1j * rng.normal(size=ns))[L]   # generic commutant element (not Hermitian)
    # identify copy 0 with each other copy through Tp (Schur: B_b^* Tp B_a is an intertwiner, a multiple of an isometry)
    ident = {}
    for dmn, js in iso.items():
        a = js[0]; ident[a] = bases_a = V[:, groups[a]]
        for b in js[1:]:
            Mab = V[:, groups[b]].conj().T @ Tp @ bases_a
            sv = np.linalg.svd(Mab, compute_uv=False)
            Mab = Mab / sv[0]
            ident[b] = V[:, groups[b]] @ Mab
    multi = {dmn: js for dmn, js in iso.items() if len(js) > 1}
    sbest = (0, None)
    if multi and k <= 16:
        for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 400):
            cols_ = []
            for dmn, js in iso.items():
                mu = len(js); r_ = rng.integers(0, mu + 1)
                if r_ == 0: continue
                if mu == 1: cols_.append(ident[js[0]]); continue
                U = rng.normal(size=(mu, r_)) + 1j * rng.normal(size=(mu, r_))
                if rng.random() < 0.3: U = np.round(U)          # rational-looking special positions too
                if np.linalg.matrix_rank(U) < r_: continue
                for c in range(r_):
                    cols_.append(sum(U[i, c] * ident[js[i]] for i in range(mu)))
            if not cols_: continue
            B = np.linalg.qr(np.hstack(cols_))[0]
            if B.shape[1] in (0, n): continue
            P = B @ B.conj().T; R = np.eye(n) - P
            fz = rk(R @ Piz @ B); fI = rk(R @ PiI @ B)
            if fz and fI / fz > sbest[0]: sbest = (fI / fz, (fI, fz, B.shape[1]))
        print("   sampled mixed submodules (multiplicities %s): best ratio %s" % ({d: len(j) for d, j in multi.items()}, sbest), flush=True)
        best = max(best, (sbest[0], None, None, None), key=lambda t: t[0])
    print("block |O|=%d |supp|=%d: stabiliser %d, commutant dim %d, constituents %d dims %s, multiplicity-free %s; "
          "best f(I)/f(zz^T) = %s (mask %s, f(I)=%s, f(move)=%s)" % (n, len(supp[Q0]), len(stab), ns, k, dims, mf,
          best[0], bin(best[1]) if best[1] else None, best[2], best[3]), flush=True)
    best_all = max(best_all, best[0])
print("m=%d: best Theorem-R ratio over submodule sums of the permutation blocks: %s (plain needs %d)" % (m, best_all, m))
