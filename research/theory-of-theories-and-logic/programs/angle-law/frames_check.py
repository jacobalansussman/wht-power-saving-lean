# Skeptic's own explicit check on the frame graph DSp(2m,2) (Lagrangians of F_2^(2m)), independent code.
# (1) build the graph by BFS under Clifford generators (symplectic action), check distance-regularity, b_j, c_j,
#     and the adjacency spectrum against theta_i = 2(2^(m-i)-1) - (2^i - 1) with the multiplicities.
# (2) Theorem F on random (non-symmetric) families of r-dim subspaces: ratio phi_j/phi_1 <= 2 - 2^(1-j).
# (3) the actual hook family Y_U = rho(g) N (N = span{e_x (x) e_T : x in T}): built by BFS over the same generators;
#     checks well-definedness, that rank excess and phi depend only on distance, and the profile vs the cell formula.
import itertools, sys
import numpy as np
from math import comb

def gens_symp(m):
    # each generator: (name, function on a 2m-bit vector (a|b) packed as (a, b) tuples)
    G = []
    for i in range(m):
        G.append(('H', i, None))
        G.append(('S', i, None))
    for c in range(m):
        for t in range(m):
            if c != t: G.append(('CX', c, t))
    return G

def act(gen, v, m):
    a, b = v
    kind, i, t = gen
    a = list(a); b = list(b)
    if kind == 'H': a[i], b[i] = b[i], a[i]
    elif kind == 'S': b[i] ^= a[i]
    else:
        c = i
        a[t] ^= a[c]; b[c] ^= b[t]
    return (tuple(a), tuple(b))

def span(vecs, m):
    S = {((0,)*m, (0,)*m)}
    for v in vecs:
        S |= {(tuple(x ^ y for x, y in zip(s[0], v[0])), tuple(x ^ y for x, y in zip(s[1], v[1]))) for s in S}
    return frozenset(S)

def lagr_act(gen, L, m):
    return frozenset(act(gen, v, m) for v in L)

def unitary(gen, m):
    n = 2**m
    U = np.zeros((n, n), dtype=complex)
    kind, i, t = gen
    for x in range(n):
        bits = [(x >> q) & 1 for q in range(m)]
        if kind == 'H':
            for y0 in (0, 1):
                bb = bits.copy(); bb[i] = y0
                y = sum(bb[q] << q for q in range(m))
                U[y, x] += (1 if not (bits[i] and y0) else -1) / np.sqrt(2)
        elif kind == 'S':
            U[x, x] = 1j if bits[i] else 1
        else:
            bb = bits.copy(); bb[t] ^= bb[i]
            y = sum(bb[q] << q for q in range(m))
            U[y, x] = 1
    return U

def wedge(g, k):
    n = g.shape[0]
    subs = list(itertools.combinations(range(n), k))
    W = np.empty((len(subs), len(subs)), dtype=complex)
    for r, T in enumerate(subs):
        gT = g[list(T), :]
        for c, S in enumerate(subs):
            W[r, c] = np.linalg.det(gT[:, list(S)])
    return W, subs

def bfs_frames(m):
    L0 = span([((0,)*m, tuple(1 if q == i else 0 for q in range(m))) for i in range(m)], m)
    gens = gens_symp(m)
    seen = {L0: 0}; order = [L0]; parent = {L0: None}
    q = [L0]
    while q:
        L = q.pop(0)
        for gi, g in enumerate(gens):
            L2 = lagr_act(g, L, m)
            if L2 not in seen:
                seen[L2] = len(order); order.append(L2); parent[L2] = (L, gi); q.append(L2)
    return order, seen, parent, gens

rng = np.random.default_rng(7)
for m in (2, 3, 4):
    order, seen, parent, gens = bfs_frames(m)
    Nf = len(order)
    dist = np.zeros((Nf, Nf), dtype=int)
    for i in range(Nf):
        for j in range(Nf):
            dist[i, j] = m - int(np.log2(len(order[i] & order[j])))
    A = (dist == 1).astype(float)
    ev = np.round(np.linalg.eigvalsh(A)).astype(int)
    vals, cnts = np.unique(ev, return_counts=True)
    th = [2*(2**(m-i)-1) - (2**i - 1) for i in range(m+1)]
    # intersection numbers from vertex 0
    ok = True
    b = []; c = []
    for j in range(m+1):
        bj = set(); cj = set()
        for v in range(Nf):
            if dist[0, v] == j:
                bj.add(int(((dist[v] == 1) & (dist[0] == j+1)).sum())); cj.add(int(((dist[v] == 1) & (dist[0] == j-1)).sum()))
        b.append(bj); c.append(cj)
    print(f"m={m}: frames {Nf}; eigenvalues {dict(zip(vals.tolist(), cnts.tolist()))}; predicted thetas {th}")
    print(f"   b_j sets {b} (pred {[2**(j+1)*(2**(m-j)-1) for j in range(m+1)]}); c_j sets {c} (pred {[2**j-1 for j in range(m+1)]})")
    # (2) random families
    if m <= 3:
        for (D, r) in [(6, 2), (10, 3), (20, 1)]:
            for trial in range(3):
                Qs = [np.linalg.qr(rng.normal(size=(D, r)) + 1j*rng.normal(size=(D, r)))[0] for _ in range(Nf)]
                Kgram = np.array([[np.linalg.norm(Qs[u].conj().T @ Qs[v])**2 for v in range(Nf)] for u in range(Nf)])
                phi = [np.mean(r - Kgram[dist == d]) for d in range(m+1)]
                rat = [phi[j]/phi[1] for j in range(1, m+1)]
                assert all(rat[j-1] <= 2 - 2**(1-j) + 1e-12 for j in range(1, m+1)), rat
            print(f"   random families D={D} r={r}: last ratios {np.round(rat, 4).tolist()} (ceilings {[2-2**(1-j) for j in range(1, m+1)]})")
    # (3) hook family
    if m <= 3:
        for k in ([1, 2] if m == 2 else [2, 4]):
            n = 2**m
            Us = [unitary(g, m) for g in gens]
            Gs = []
            for U in Us:
                W, subs = wedge(U, k); Gs.append(np.kron(U, W))
            inN = np.array([x in T for x in range(n) for T in subs])
            Q0 = np.eye(len(inN), dtype=complex)[:, inN]
            Q = {order[0]: Q0}
            for L in order[1:]:
                Lp, gi = parent[L]
                Q[L] = Gs[gi] @ Q[Lp]
            # well-definedness: for every frame and generator, G Q[L] spans the same space as Q[gL]
            worst = 0
            for L in order:
                for gi, g in enumerate(gens):
                    L2 = lagr_act(g, L, m)
                    A1 = Gs[gi] @ Q[L]; B1 = Q[L2]
                    worst = max(worst, np.linalg.norm(A1 - B1 @ (B1.conj().T @ A1)))
            rdim = Q0.shape[1]
            rk = {d: set() for d in range(m+1)}; ph = {d: [] for d in range(m+1)}
            pairs = [(0, v) for v in range(Nf)] + [(int(rng.integers(Nf)), int(rng.integers(Nf))) for _ in range(40)]
            for (u, v) in pairs:
                d = dist[u, v]
                Qu, Qv = Q[order[u]], Q[order[v]]
                Mx = Qv - Qu @ (Qu.conj().T @ Qv)
                s = np.linalg.svd(Mx, compute_uv=False)
                rk[d].add(int((s > 1e-8).sum())); ph[d].append(rdim - np.linalg.norm(Qu.conj().T @ Qv)**2)
            phs = {d: (round(min(ph[d]), 6), round(max(ph[d]), 6)) for d in ph}
            print(f"   hook k={k}: dim N {rdim}; well-defined (max residual {worst:.1e}); rank excess by distance {rk}; phi (min,max) {phs}")
            p1 = np.mean(ph[1])
            print(f"     phi ratios {[round(np.mean(ph[j])/p1, 6) for j in range(1, m+1)]}; closed form phi_j = {[comb(n-2,k-1)*n*(1-2**-j) for j in range(1, m+1)]}")
