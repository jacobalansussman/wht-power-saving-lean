# nclp.py -- entropy-LP lower bounds on the cost of a network code (min-cost form), with orbit variables.
# Network: nodes, directed arcs (an undirected edge = two arcs, paid separately), sessions (source node, sink nodes,
# rate). A code: every node stores everything it has seen and computes freely; each use of an arc carries one symbol.
# Cost = arc uses >= sum over arcs of h(X_arc) (X_arc = everything ever sent on the arc).
# Ground set: one element per session (Y_i) and per arc (X_a), or per GROUP of them (a group is the union; valid).
# Constraints (all valid for every code; linear rank functions satisfy all of them):
#   polymatroid (elemental Shannon inequalities), h(empty) = 0;
#   sources independent with h(Y_i) = rate_i  (or rates free with a normalisation, for the worst case);
#   everything sent is a function of all sources;
#   for every node set Q: h(In(Q)) = h(In(Q) + Out(Q) + Dec(Q)), In(Q) = sources at Q + arcs entering Q,
#   Out(Q) = arcs leaving Q, Dec(Q) = sessions with a sink in Q  (causality: all Q knows came in or started there).
# Optional cutting planes: Ingleton (all fields) and DFZ inequality (65) of arXiv:1311.4601 Thm 8.6 (char != 2).
import itertools, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

class Net:
    def __init__(self, nodes, arcs, sessions):
        """nodes: list of names; arcs: list of (u, v) (each a separate ground element); sessions: list of
        (src, [sinks], rate or None)."""
        self.nodes = list(nodes); self.arcs = list(arcs); self.sessions = list(sessions)
        self.ns = len(self.sessions); self.na = len(self.arcs); self.n = self.ns + self.na
    def bit_session(self, i): return 1 << i
    def bit_arc(self, a): return 1 << (self.ns + a)

def orbit_map(n, perms):
    """perms: list of permutations of range(n) (a generating set is NOT enough: pass the whole group).
    returns (canon: array mapping mask -> orbit id, nvars)"""
    N = 1 << n
    masks = np.arange(N, dtype=np.int64)
    best = masks.copy()
    for p in perms:
        img = np.zeros(N, dtype=np.int64)
        for i in range(n):
            img |= ((masks >> i) & 1) << p[i]
        best = np.minimum(best, img)
    reps, canon = np.unique(best, return_inverse=True)
    return canon.astype(np.int64), len(reps)

def close_group(gens, n):
    G = {tuple(range(n))}; frontier = [tuple(range(n))]
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = tuple(s[g[i]] for i in range(n))
                if h not in G: G.add(h); new.append(h)
        frontier = new
    return sorted(G)

def build_lp(net, perms=None, rates_free=False, routing=None, extra_rows=None, verbose=True):
    n = net.n; N = 1 << n
    t0 = time.time()
    if perms:
        canon, nv = orbit_map(n, perms)
    else:
        canon, nv = np.arange(N, dtype=np.int64), N
    full = N - 1
    rows_ub = []   # (cols array, vals array), meaning sum vals*h[cols] <= 0  (homogeneous) -- built vectorised
    # elemental: monotone h(N - i) <= h(N)
    I = []; J = []; V = []; r = 0
    for i in range(n):
        I += [r, r]; J += [canon[full & ~(1 << i)], canon[full]]; V += [1.0, -1.0]; r += 1
    # submodularity: h(S+i+j) + h(S) - h(S+i) - h(S+j) <= 0, S subset of rest
    masks = np.arange(N, dtype=np.int64)
    sub_blocks = []
    for i in range(n):
        for j in range(i + 1, n):
            S = masks[((masks >> i) & 1) == 0]
            S = S[((S >> j) & 1) == 0]
            a = canon[S | (1 << i) | (1 << j)]; b = canon[S]; c = canon[S | (1 << i)]; d = canon[S | (1 << j)]
            sub_blocks.append(np.stack([a, b, c, d], axis=1))
    sub = np.concatenate(sub_blocks, axis=0)
    # canonicalise: (a,b) as an unordered pair with +1, (c,d) unordered with -1; drop rows that are 0
    p1 = np.sort(sub[:, :2], axis=1); p2 = np.sort(sub[:, 2:], axis=1)
    key = np.concatenate([p1, p2], axis=1)
    key = np.unique(key, axis=0)
    # remove trivial rows where {a,b} == {c,d}
    nontriv = ~((key[:, 0] == key[:, 2]) & (key[:, 1] == key[:, 3]))
    key = key[nontriv]
    m_sub = key.shape[0]
    rr = np.arange(r, r + m_sub)
    I = np.concatenate([np.array(I, dtype=np.int64), rr, rr, rr, rr])
    J = np.concatenate([np.array(J, dtype=np.int64), key[:, 0], key[:, 1], key[:, 2], key[:, 3]])
    V = np.concatenate([np.array(V), np.ones(m_sub), np.ones(m_sub), -np.ones(m_sub), -np.ones(m_sub)])
    r += m_sub
    n_ub = r
    A_ub = sp.csr_matrix((V, (I, J)), shape=(n_ub, nv))   # duplicates summed (e.g. a == b)
    b_ub = np.zeros(n_ub)
    # equalities
    E = []   # list of dict col->coef, rhs
    def eq(d, rhs): E.append((d, rhs))
    def add(d, mask, c):
        k = int(canon[mask]); d[k] = d.get(k, 0.0) + c
    eq({int(canon[0]): 1.0}, 0.0)
    smask = (1 << net.ns) - 1
    # sources independent: h(Y_S) = sum h(Y_i); rates
    for k in range(1, net.ns + 1):
        for S in itertools.combinations(range(net.ns), k):
            if k == 1: continue
            d = {}; add(d, sum(1 << i for i in S), 1.0)
            for i in S: add(d, 1 << i, -1.0)
            eq(d, 0.0)
    if not rates_free:
        for i, (src, sinks, rate) in enumerate(net.sessions):
            eq({int(canon[1 << i]): 1.0}, float(rate))
    else:
        # normalise routing cost: sum_i routing_i * h(Y_i) = 1
        d = {}
        for i in range(net.ns): add(d, 1 << i, float(routing[i]))
        eq(d, 1.0)
    # everything is a function of the sources
    d = {}; add(d, full, 1.0); add(d, smask, -1.0); eq(d, 0.0)
    # closure for every node set Q
    nodes = net.nodes; nn = len(nodes)
    for q in range(1, 1 << nn):
        Q = {nodes[t] for t in range(nn) if q >> t & 1}
        inm = 0; full_m = 0
        for i, (src, sinks, rate) in enumerate(net.sessions):
            if src in Q: inm |= 1 << i
            if any(s in Q for s in sinks): full_m |= 1 << i
        full_m |= inm
        for a, (u, v) in enumerate(net.arcs):
            if (u not in Q) and (v in Q): inm |= net.bit_arc(a)
            if (u in Q) and (v not in Q): full_m |= net.bit_arc(a)
        full_m |= inm
        if full_m != inm:
            d = {}; add(d, full_m, 1.0); add(d, inm, -1.0)
            if d and any(abs(x) > 0 for x in d.values()): eq(d, 0.0)
    if extra_rows:
        for (d0, rhs) in extra_rows(canon): eq(d0, rhs)
    # dedupe equalities
    seen = set(); E2 = []
    for d, rhs in E:
        d = {k: v for k, v in d.items() if abs(v) > 1e-15}
        if not d: continue
        key2 = (tuple(sorted(d.items())), rhs)
        if key2 in seen: continue
        seen.add(key2); E2.append((d, rhs))
    Ie = []; Je = []; Ve = []; be = []
    for t, (d, rhs) in enumerate(E2):
        for k, v in d.items(): Ie.append(t); Je.append(k); Ve.append(v)
        be.append(rhs)
    A_eq = sp.csr_matrix((Ve, (Ie, Je)), shape=(len(E2), nv)); b_eq = np.array(be)
    # objective: sum over arcs h(X_a)
    c = np.zeros(nv)
    for a in range(net.na): c[canon[net.bit_arc(a)]] += 1.0
    if verbose:
        print(f"  LP: ground {n}, variables {nv}, inequalities {n_ub}, equalities {len(E2)}, built in {time.time()-t0:.1f}s")
        sys.stdout.flush()
    return dict(c=c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, canon=canon, nv=nv, n=n)

def solve(lp, extra_ub=None, verbose=True, method='highs'):
    A_ub, b_ub = lp['A_ub'], lp['b_ub']
    if extra_ub:
        rows = []; 
        for d in extra_ub:
            rows.append(d)
        Ix = []; Jx = []; Vx = []
        for t, d in enumerate(rows):
            for k, v in d.items(): Ix.append(t); Jx.append(k); Vx.append(v)
        X = sp.csr_matrix((Vx, (Ix, Jx)), shape=(len(rows), lp['nv']))
        A_ub = sp.vstack([A_ub, X]).tocsr(); b_ub = np.concatenate([b_ub, np.zeros(len(rows))])
    t0 = time.time()
    res = linprog(lp['c'], A_ub=A_ub, b_ub=b_ub, A_eq=lp['A_eq'], b_eq=lp['b_eq'], bounds=(0, None), method=method)
    if verbose:
        print(f"  solved: status {res.status} ({res.message[:40]}), value {res.fun if res.status==0 else None}, {time.time()-t0:.1f}s")
        sys.stdout.flush()
    return res

# ---- linear rank inequalities as cutting planes (letters = arbitrary unions of ground elements) ----
def H(h, canon, mask): return h[canon[mask]]

def ingleton_row(canon, A, B, C, D):
    """Ingleton: h(A)+h(B)+h(ACD... ) form:  I(A;B) <= I(A;B|C) + I(A;B|D) + I(C;D), written as  lhs - rhs <= 0
    I(A;B) = hA + hB - hAB ; I(A;B|C) = hAC + hBC - hABC - hC ; I(C;D) = hC + hD - hCD."""
    d = {}
    def add(m, v):
        k = int(canon[m]); d[k] = d.get(k, 0.0) + v
    add(A, 1); add(B, 1); add(A | B, -1)
    add(A | C, -1); add(B | C, -1); add(A | B | C, 1); add(C, 1)
    add(A | D, -1); add(B | D, -1); add(A | B | D, 1); add(D, 1)
    add(C, -1); add(D, -1); add(C | D, 1)
    return {k: v for k, v in d.items() if abs(v) > 1e-12}

def dfz65_row(canon, L):
    """DFZ (65), arXiv:1311.4601 Thm 8.6, valid for linear rank functions over fields of characteristic != 2:
    2H(A)+H(B)+2H(C) <= H(W)+H(X)+H(Y)+H(Z)+2H(A|Z,Y)+H(B|X,Z)+2H(C|A,X)+3H(X|W,Y)+3H(Z|W,C)+5H(W|A,B)
                        +5H(Y|B,C)+5(H(A)+H(B)+H(C)-H(A,B,C)).   Returned as lhs - rhs <= 0."""
    A, B, C, W, X, Y, Z = L
    d = {}
    def add(m, v):
        k = int(canon[m]); d[k] = d.get(k, 0.0) + v
    def Hc(x, ys, coef):  # coef * H(x | ys) = coef*(H(x,ys) - H(ys))
        yy = 0
        for y in ys: yy |= y
        add(x | yy, coef); add(yy, -coef)
    # lhs
    add(A, 2); add(B, 1); add(C, 2)
    # minus rhs
    for v in (W, X, Y, Z): add(v, -1)
    Hc(A, [Z, Y], -2); Hc(B, [X, Z], -1); Hc(C, [A, X], -2); Hc(X, [W, Y], -3); Hc(Z, [W, C], -3)
    Hc(W, [A, B], -5); Hc(Y, [B, C], -5)
    add(A, -5); add(B, -5); add(C, -5); add(A | B | C, 5)
    return {k: v for k, v in d.items() if abs(v) > 1e-12}

def row_value(d, h): return sum(v * h[k] for k, v in d.items())

def find_cuts(h, canon, n, letters, kind, max_new=200, tries=20000, rng=None, tol=1e-7):
    """search violated Ingleton / DFZ instances with letters drawn from `letters` (list of masks)."""
    import random
    rng = rng or random.Random(0)
    found = {}
    if kind == 'ingleton':
        # exhaustive over 4-tuples of letters if small, else random
        L = letters
        total = len(L) ** 4
        it = itertools.product(range(len(L)), repeat=4) if total <= tries else (tuple(rng.randrange(len(L)) for _ in range(4)) for _ in range(tries))
        for (a, b, c, e) in it:
            A, B, C, D = L[a], L[b], L[c], L[e]
            if len({A, B, C, D}) < 4: continue
            d = ingleton_row(canon, A, B, C, D)
            v = row_value(d, h)
            if v > tol:
                key = tuple(sorted(d.items()))
                if key not in found: found[key] = (v, d)
    else:
        for _ in range(tries):
            Ls = [letters[rng.randrange(len(letters))] for _ in range(7)]
            if len(set(Ls)) < 7: continue
            d = dfz65_row(canon, Ls)
            v = row_value(d, h)
            if v > tol:
                key = tuple(sorted(d.items()))
                if key not in found: found[key] = (v, d)
    best = sorted(found.values(), key=lambda t: -t[0])[:max_new]
    return [d for v, d in best], (best[0][0] if best else 0.0)

# ---- cutting-plane version: same LP, elemental submodularity rows added lazily (exact same optimum) ----
def elemental_arrays(n, canon):
    masks = np.arange(1 << n, dtype=np.int64)
    blocks = []
    for i in range(n):
        for j in range(i + 1, n):
            S = masks[(((masks >> i) & 1) == 0) & (((masks >> j) & 1) == 0)]
            blocks.append(np.stack([canon[S | (1 << i) | (1 << j)], canon[S], canon[S | (1 << i)], canon[S | (1 << j)]], axis=1))
    sub = np.concatenate(blocks, axis=0)
    p1 = np.sort(sub[:, :2], axis=1); p2 = np.sort(sub[:, 2:], axis=1)
    key = np.unique(np.concatenate([p1, p2], axis=1), axis=0)
    key = key[~((key[:, 0] == key[:, 2]) & (key[:, 1] == key[:, 3]))]
    return key

def solve_cp(net, perms=None, rates_free=False, routing=None, batch=40000, tol=1e-9, verbose=True, method='highs',
             extra_ub=None, max_rounds=200):
    """cutting-plane solve of the full Shannon LP: start from monotonicity + network rows, add violated elemental
    submodularity rows in batches until none is violated (then the optimum equals the full LP's)."""
    t0 = time.time()
    n = net.n
    if perms: canon, nv = orbit_map(n, perms)
    else: canon, nv = np.arange(1 << n, dtype=np.int64), 1 << n
    key = elemental_arrays(n, canon)
    A_eq, b_eq, c = _eq_part(net, canon, nv, rates_free, routing)
    full = (1 << n) - 1
    mono = []
    for i in range(n):
        mono.append({int(canon[full & ~(1 << i)]): 1.0, int(canon[full]): -1.0})
    active = np.zeros(len(key), dtype=bool)
    # seed: rows with S = empty set (pairwise submodularity of singletons) are cheap and useful
    rounds = 0
    while True:
        rows = mono[:] + ([] if extra_ub is None else list(extra_ub))
        K = key[active]
        I = np.repeat(np.arange(len(K)), 4); J = K.reshape(-1); V = np.tile(np.array([1.0, 1.0, -1.0, -1.0]), len(K))
        A1 = sp.csr_matrix((V, (I, J)), shape=(len(K), nv))
        Im = []; Jm = []; Vm = []
        for t, d in enumerate(rows):
            for k, v in d.items(): Im.append(t); Jm.append(k); Vm.append(v)
        A0 = sp.csr_matrix((Vm, (Im, Jm)), shape=(len(rows), nv))
        A_ub = sp.vstack([A0, A1]).tocsr()
        res = linprog(c, A_ub=A_ub, b_ub=np.zeros(A_ub.shape[0]), A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method=method)
        if res.status != 0:
            print("  LP status", res.status, res.message); return res, canon, key, active
        h = res.x
        viol = h[key[:, 0]] + h[key[:, 1]] - h[key[:, 2]] - h[key[:, 3]]
        bad = np.nonzero((viol > tol) & ~active)[0]
        rounds += 1
        if verbose:
            print(f"  cp round {rounds}: value {res.fun:.9f}, active rows {active.sum()}, violated {len(bad)} (max {viol.max():.2e}), {time.time()-t0:.0f}s")
            sys.stdout.flush()
        if len(bad) == 0 or rounds >= max_rounds:
            res.canon = canon
            return res, canon, key, active
        order = bad[np.argsort(-viol[bad])][:batch]
        active[order] = True

def _eq_part(net, canon, nv, rates_free, routing):
    n = net.n; full = (1 << n) - 1
    E = []
    def eq(d, rhs): E.append((d, rhs))
    def add(d, mask, c0):
        k = int(canon[mask]); d[k] = d.get(k, 0.0) + c0
    eq({int(canon[0]): 1.0}, 0.0)
    smask = (1 << net.ns) - 1
    for k in range(2, net.ns + 1):
        for S in itertools.combinations(range(net.ns), k):
            d = {}; add(d, sum(1 << i for i in S), 1.0)
            for i in S: add(d, 1 << i, -1.0)
            eq(d, 0.0)
    if not rates_free:
        for i, (src, sinks, rate) in enumerate(net.sessions): eq({int(canon[1 << i]): 1.0}, float(rate))
    else:
        d = {}
        for i in range(net.ns): add(d, 1 << i, float(routing[i]))
        eq(d, 1.0)
    d = {}; add(d, full, 1.0); add(d, smask, -1.0); eq(d, 0.0)
    nodes = net.nodes; nn = len(nodes)
    for q in range(1, 1 << nn):
        Q = {nodes[t] for t in range(nn) if q >> t & 1}
        inm = 0; fm = 0
        for i, (src, sinks, rate) in enumerate(net.sessions):
            if src in Q: inm |= 1 << i
            if any(s in Q for s in sinks): fm |= 1 << i
        fm |= inm
        for a, (u, v) in enumerate(net.arcs):
            if (u not in Q) and (v in Q): inm |= net.bit_arc(a)
            if (u in Q) and (v not in Q): fm |= net.bit_arc(a)
        fm |= inm
        if fm != inm:
            d = {}; add(d, fm, 1.0); add(d, inm, -1.0); eq(d, 0.0)
    seen = set(); E2 = []
    for d, rhs in E:
        d = {k: v for k, v in d.items() if abs(v) > 1e-15}
        if not d: continue
        kk = (tuple(sorted(d.items())), rhs)
        if kk in seen: continue
        seen.add(kk); E2.append((d, rhs))
    Ie = []; Je = []; Ve = []; be = []
    for t, (d, rhs) in enumerate(E2):
        for k, v in d.items(): Ie.append(t); Je.append(k); Ve.append(v)
        be.append(rhs)
    A_eq = sp.csr_matrix((Ve, (Ie, Je)), shape=(len(E2), nv)); b_eq = np.array(be)
    c = np.zeros(nv)
    for a in range(net.na): c[canon[net.bit_arc(a)]] += 1.0
    return A_eq, b_eq, c

# ---- closure-reduced LP: identify every set with its functional closure (valid: h(S) = h(cl S) for every code) ----
def closure_rules(net):
    """list of (premise mask, conclusion mask): if premise is inside S then conclusion joins S."""
    rules = []
    nodes = net.nodes; nn = len(nodes)
    for q in range(1, 1 << nn):
        Q = {nodes[t] for t in range(nn) if q >> t & 1}
        inm = 0; fm = 0
        for i, (src, sinks, rate) in enumerate(net.sessions):
            if src in Q: inm |= 1 << i
            if any(s in Q for s in sinks): fm |= 1 << i
        for a, (u, v) in enumerate(net.arcs):
            if (u not in Q) and (v in Q): inm |= net.bit_arc(a)
            if (u in Q) and (v not in Q): fm |= net.bit_arc(a)
        fm |= inm
        if fm != inm: rules.append((inm, fm))
    smask = (1 << net.ns) - 1; full = (1 << net.n) - 1
    rules.append((smask, full))
    return rules

def closure_map(net):
    N = 1 << net.n
    cl = np.arange(N, dtype=np.int64)
    rules = closure_rules(net)
    changed = True
    while changed:
        changed = False
        for (p, c) in rules:
            hit = (cl & p) == p
            new = np.where(hit, cl | c, cl)
            if (new != cl).any(): changed = True; cl = new
    return cl

def build_closed_lp(net, perms=None, rates_free=False, routing=None, verbose=True):
    n = net.n; N = 1 << n
    t0 = time.time()
    cl = closure_map(net)
    if perms: canon0, _ = orbit_map(n, perms)
    else: canon0 = np.arange(N, dtype=np.int64)
    key0 = canon0[cl]                      # orbit id of the closure
    reps, canon = np.unique(key0, return_inverse=True); canon = canon.astype(np.int64); nv = len(reps)
    full = N - 1
    masks = np.arange(N, dtype=np.int64)
    blocks = []
    for i in range(n):
        for j in range(i + 1, n):
            S = masks[(((masks >> i) & 1) == 0) & (((masks >> j) & 1) == 0)]
            blocks.append(np.stack([canon[S | (1 << i) | (1 << j)], canon[S], canon[S | (1 << i)], canon[S | (1 << j)]], axis=1))
    sub = np.concatenate(blocks, axis=0)
    p1 = np.sort(sub[:, :2], axis=1); p2 = np.sort(sub[:, 2:], axis=1)
    key = np.unique(np.concatenate([p1, p2], axis=1), axis=0)
    key = key[~((key[:, 0] == key[:, 2]) & (key[:, 1] == key[:, 3]))]
    mono = np.unique(np.stack([canon[full & ~(1 << i)] for i in range(n)]))
    m_sub = key.shape[0]
    rr = np.arange(m_sub)
    I = np.concatenate([rr, rr, rr, rr]); J = np.concatenate([key[:, 0], key[:, 1], key[:, 2], key[:, 3]])
    V = np.concatenate([np.ones(m_sub), np.ones(m_sub), -np.ones(m_sub), -np.ones(m_sub)])
    # monotonicity h(N - i) <= h(N)
    r0 = m_sub; Im = []; Jm = []; Vm = []
    for t, k in enumerate(mono):
        if k != canon[full]:
            Im += [r0, r0]; Jm += [int(k), int(canon[full])]; Vm += [1.0, -1.0]; r0 += 1
    I = np.concatenate([I, np.array(Im, dtype=np.int64)]); J = np.concatenate([J, np.array(Jm, dtype=np.int64)])
    V = np.concatenate([V, np.array(Vm)])
    A_ub = sp.csr_matrix((V, (I, J)), shape=(r0, nv))
    # equalities: empty set, source independence, rates (closure rows are now identities)
    E = []
    def add(d, mask, c0):
        k = int(canon[mask]); d[k] = d.get(k, 0.0) + c0
    E.append(({int(canon[0]): 1.0}, 0.0))
    for k in range(2, net.ns + 1):
        for S in itertools.combinations(range(net.ns), k):
            d = {}; add(d, sum(1 << i for i in S), 1.0)
            for i in S: add(d, 1 << i, -1.0)
            E.append((d, 0.0))
    if not rates_free:
        for i, (src, sinks, rate) in enumerate(net.sessions): E.append(({int(canon[1 << i]): 1.0}, float(rate)))
    else:
        d = {}
        for i in range(net.ns): add(d, 1 << i, float(routing[i]))
        E.append((d, 1.0))
    seen = set(); E2 = []
    for d, rhs in E:
        d = {k: v for k, v in d.items() if abs(v) > 1e-15}
        if not d:
            if abs(rhs) > 1e-12: raise ValueError("inconsistent equality after closure")
            continue
        kk = (tuple(sorted(d.items())), rhs)
        if kk not in seen: seen.add(kk); E2.append((d, rhs))
    Ie = []; Je = []; Ve = []; be = []
    for t, (d, rhs) in enumerate(E2):
        for k, v in d.items(): Ie.append(t); Je.append(k); Ve.append(v)
        be.append(rhs)
    A_eq = sp.csr_matrix((Ve, (Ie, Je)), shape=(len(E2), nv)); b_eq = np.array(be)
    c = np.zeros(nv)
    for a in range(net.na): c[canon[net.bit_arc(a)]] += 1.0
    if verbose:
        print(f"  closed LP: ground {n}, closed-set orbit variables {nv} (of {N}), inequalities {A_ub.shape[0]}, equalities {len(E2)}, {time.time()-t0:.1f}s")
        sys.stdout.flush()
    return dict(c=c, A_ub=A_ub, b_ub=np.zeros(A_ub.shape[0]), A_eq=A_eq, b_eq=b_eq, canon=canon, nv=nv, n=n)

# ---- first-order solve (HiGHS PDLP via highspy) with a RIGOROUS lower bound from the approximate dual ----
def solve_pdlp(lp, Rtot, extra_ub=None, time_limit=3000, tol=1e-8, verbose=True):
    """min c.x, A_ub x <= 0, A_eq x = b_eq, x >= 0. Every variable is some h(S) <= h(all) = Rtot (total source rate),
    so for ANY row multipliers y (ub rows clipped to <= 0) the number
        LB = b_eq . y_eq + sum_j min(0, (c - A^T y)_j) * Rtot
    is a valid lower bound on the LP optimum (weak duality with box 0 <= x <= Rtot). Returns (primal, LB, res)."""
    import highspy
    A_ub, b_ub = lp['A_ub'], lp['b_ub']
    if extra_ub:
        Ix = []; Jx = []; Vx = []
        for t, d in enumerate(extra_ub):
            for k, v in d.items(): Ix.append(t); Jx.append(k); Vx.append(v)
        X = sp.csr_matrix((Vx, (Ix, Jx)), shape=(len(extra_ub), lp['nv']))
        A_ub = sp.vstack([A_ub, X]).tocsr(); b_ub = np.concatenate([b_ub, np.zeros(len(extra_ub))])
    A = sp.vstack([A_ub, lp['A_eq']]).tocsc()
    nu = A_ub.shape[0]
    lo = np.concatenate([np.full(nu, -highspy.kHighsInf), lp['b_eq']]); hi = np.concatenate([b_ub, lp['b_eq']])
    h = highspy.Highs()
    h.setOptionValue('output_flag', False); h.setOptionValue('solver', 'pdlp'); h.setOptionValue('time_limit', float(time_limit))
    h.setOptionValue('threads', 1)
    h.setOptionValue('primal_feasibility_tolerance', tol); h.setOptionValue('dual_feasibility_tolerance', tol)
    model = highspy.HighsLp()
    model.num_col_ = lp['nv']; model.num_row_ = A.shape[0]
    model.col_cost_ = lp['c']; model.col_lower_ = np.zeros(lp['nv']); model.col_upper_ = np.full(lp['nv'], float(Rtot))
    model.row_lower_ = lo; model.row_upper_ = hi
    model.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    model.a_matrix_.start_ = A.indptr; model.a_matrix_.index_ = A.indices; model.a_matrix_.value_ = A.data
    h.passModel(model)
    t0 = time.time(); h.run()
    sol = h.getSolution(); info = h.getInfo()
    y = np.array(sol.row_dual); x = np.array(sol.col_value)
    y_ub = np.minimum(y[:nu], 0.0); y_eq = y[nu:]
    red = lp['c'] - A_ub.T @ y_ub - lp['A_eq'].T @ y_eq
    LB = float(lp['b_eq'] @ y_eq + np.minimum(red, 0.0).sum() * Rtot)
    primal = float(lp['c'] @ x)
    if verbose:
        print(f"  pdlp: status {h.modelStatusToString(h.getModelStatus())}, primal {primal:.9f}, rigorous dual bound {LB:.9f}, {time.time()-t0:.0f}s")
        sys.stdout.flush()
    class R: pass
    r = R(); r.x = x; r.fun = primal; r.lb = LB; r.y = y
    return primal, LB, r
