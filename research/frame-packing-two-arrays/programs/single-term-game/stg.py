#!/usr/bin/env python3
"""stg.py (p6-refute, session 2, 15 Oct 2026): the SINGLE-TERM two-input module game, exact and exhaustive.

Class (p6-prove PROOF.md section 4): two inputs, any number of scratch registers, every MOVED register holds at most
one term per input.  By Lemma M + Lemma G this is a game on positions (frames) only:
  - a single of input l at position x (available position);
  - a pair = lineage copy at input-1 position x, input-2 position tau[x] (tau a twist = right multiplication by a
    group element, acting on the frames); moving a pair moves x to any neighbour x2 and y to tau[x2];
  - availability: a component of the bipartite position graph (vertices (1,x), (2,y); one edge per locked lineage
    copy) is UNLOCKED if it contains an available vertex or a cycle (permissive closure, any cycle unlocks: an
    over-approximation of Lemma G's "non-degenerate cycle", the safe direction for lower bounds);
  - moves (cost 1): a single; a lineage copy (same lineage); a CREATED pair (both ends available, any twist, new
    lineage); a COMPOSITE (path of >= 3 edges in a locked tree component, twist = c o b^-1 o a, new lineage).
Checked at every state: (I) max level(A1) + max level(A2) <= M  (Conjecture T), and the goal (both inputs available at
a target-orbit frame).  Symmetry: independent right multiplication by free elements on each input, and the swap.

Games: any toylab FrameGame (positions = right cosets P g, moves P g -> P h g, twists P g -> P g t), or Sp(2m, 2)
built from p6-prove's frames.py (positions = Lagrangians).
"""
import sys, time, itertools
from collections import defaultdict

class Game:
    def __init__(self, N, nb, lv, targets, twists, free_idx, name="", mult_table=True):
        """N positions; nb[x] neighbour lists; lv[x] level (distance from the raw frame, which has lv 0);
        targets: set of target positions (a free orbit at level m); twists: list of permutations (tuples), index 0 the
        identity, closed under composition and inverse; free_idx: indices of the twists fixing the raw frame."""
        self.N, self.nb, self.lv, self.name = N, nb, lv, name
        self.targets = set(targets)
        self.tw = twists; self.ti = {t: i for i, t in enumerate(twists)}
        assert len(self.ti) == len(twists)
        self.T = len(twists)
        self.inv = [self.ti[tuple(sorted(range(N), key=lambda i: t[i]))] for t in twists]
        self.free = free_idx
        self.m = max(lv[x] for x in targets)
        self.tab = None
        if mult_table:
            self.tab = [[self.ti[tuple(a[b[i]] for i in range(N))] for b in twists] for a in twists]
        # for created pairs: by (x0) the twists with tau[x0] = y, grouped by y
        self.by_x_y = {}
        self.maxlv = max(lv)

    def comp(self, a, b):                       # twist a after b
        if self.tab is not None: return self.tab[a][b]
        A, B = self.tw[a], self.tw[b]
        return self.ti[tuple(A[B[i]] for i in range(self.N))]

    def twists_from(self, x0, y):
        k = (x0, y)
        if k not in self.by_x_y:
            self.by_x_y[k] = [i for i, t in enumerate(self.tw) if t[x0] == y]
        return self.by_x_y[k]

# ---------------------------------------------------------------------------------------------- state machinery
def bits(S):
    while S:
        b = S & -S; yield b.bit_length() - 1; S ^= b

def maxlv(g, S):
    return max((g.lv[x] for x in bits(S)), default=-1)

def closure(g, A1, A2, lins):
    """permissive closure; returns (A1, A2, lins) with every copy lying in an unlocked component removed"""
    N = g.N
    while True:
        par = {}
        def f(a):
            while par.get(a, a) != a:
                par[a] = par.get(par[a], par[a]); a = par[a]
            return a
        edges = []
        for (t, xs) in lins:
            tt = g.tw[t]
            for x in bits(xs):
                u, w = x, N + tt[x]
                edges.append((u, w))
                ru, rw = f(u), f(w)
                if ru != rw: par[ru] = rw
        if not edges: return A1, A2, lins
        V = defaultdict(set); E = defaultdict(int)
        for (u, w) in edges:
            r = f(u); V[r].add(u); V[r].add(w); E[r] += 1
        nA1, nA2 = A1, A2
        for r, vs in V.items():
            avail = any(((A1 >> v) & 1) if v < N else ((A2 >> (v - N)) & 1) for v in vs)
            if avail or E[r] >= len(vs):
                for v in vs:
                    if v < N: nA1 |= 1 << v
                    else: nA2 |= 1 << (v - N)
        nl = []
        for (t, xs) in lins:
            xs2 = xs & ~nA1
            if xs2: nl.append((t, xs2))
        nl = tuple(sorted(nl))
        if nA1 == A1 and nA2 == A2 and nl == lins: return A1, A2, lins
        A1, A2, lins = nA1, nA2, nl

def apply_mask(p, S):
    r = 0
    for x in bits(S): r |= 1 << p[x]
    return r

def canon(g, A1, A2, lins):
    best = None
    for sw in (0, 1):
        if sw:
            B1, B2 = A2, A1
            L = []
            for (t, xs) in lins:
                tt = g.tw[t]
                L.append((g.inv[t], apply_mask(tt, xs)))
        else:
            B1, B2, L = A1, A2, list(lins)
        # argmin over free a of a(B1); over free b of b(B2)
        ims1 = [(apply_mask(g.tw[a], B1), a) for a in g.free]
        m1 = min(v for v, _ in ims1); As = [a for v, a in ims1 if v == m1]
        ims2 = [(apply_mask(g.tw[b], B2), b) for b in g.free]
        m2 = min(v for v, _ in ims2); Bs = [b for v, b in ims2 if v == m2]
        if best is not None and (m1, m2) > best[:2]: continue
        for a in As:
            ai = g.inv[a]; pa = g.tw[a]
            for b in Bs:
                LL = tuple(sorted((g.comp(b, g.comp(t, ai)), apply_mask(pa, xs)) for (t, xs) in L))
                key = (m1, m2, LL)
                if best is None or key < best: best = key
    return best

def components_tree_paths(g, A1, A2, lins):
    """composites: for every locked (tree) component, every side-1 vertex u and side-2 vertex w joined by a path of
    >= 3 edges: (twist u -> w, u)"""
    N = g.N
    adj = defaultdict(list)
    for (t, xs) in lins:
        tt = g.tw[t]
        for x in bits(xs):
            u, w = x, N + tt[x]
            adj[u].append((w, t)); adj[w].append((u, t))
    out = set()
    for start in [v for v in adj if v < N]:
        # DFS: twist mapping start position -> current vertex position (side 2), or side-1 -> side-1 map
        tw = {start: (0, 0)}       # vertex -> (twist index, number of edges)
        st = [start]
        while st:
            v = st.pop()
            t, k = tw[v]
            for (nbv, s) in adj[v]:
                if nbv in tw: continue
                if v < N: nt = g.comp(s, t)            # side 1 -> side 2 along edge with twist s
                else: nt = g.comp(g.inv[s], t)         # side 2 -> side 1
                tw[nbv] = (nt, k + 1); st.append(nbv)
                if nbv >= N and k + 1 >= 3: out.add((nt, start))
    return out

def successors(g, st, with_ends=False):
    A1, A2, lins = st
    N = g.N
    res = []
    for l in (1, 2):
        A = A1 if l == 1 else A2
        for x in bits(A):
            for x2 in g.nb[x]:
                if not (A >> x2) & 1:
                    res.append(("s%d" % l, (A1 | (1 << x2), A2, lins) if l == 1 else (A1, A2 | (1 << x2), lins), [(l - 1, x2)]))
    for i, (t, xs) in enumerate(lins):
        for x in bits(xs):
            for x2 in g.nb[x]:
                if not (xs >> x2) & 1:
                    nl = list(lins); nl[i] = (t, xs | (1 << x2))
                    res.append(("c", (A1, A2, tuple(sorted(nl))), [(0, x2), (1, g.tw[t][x2])]))
    for x0 in bits(A1):
        for y0 in bits(A2):
            for t in g.twists_from(x0, y0):
                tt = g.tw[t]
                for x2 in g.nb[x0]:
                    if (A1 >> x2) & 1 or (A2 >> tt[x2]) & 1: continue      # dominated by a single move
                    res.append(("p", (A1, A2, tuple(sorted(lins + ((t, 1 << x2),)))), [(0, x2), (1, tt[x2])]))
    for (t, u) in components_tree_paths(g, A1, A2, lins):
        for u2 in g.nb[u]:
            res.append(("k", (A1, A2, tuple(sorted(lins + ((t, 1 << u2),)))), [(0, u2), (1, g.tw[t][u2])]))
    return res if with_ends else [(a, b) for (a, b, c) in res]

def bfs(g, D, report=print, keep_examples=True, max_states=None):
    A1 = A2 = 1 << 0                         # raw frame is position 0
    st = closure(g, A1, A2, ())
    layer = {canon(g, *st): st}
    parent = {}
    hist = []
    for M in range(1, D + 1):
        t0 = time.time()
        new = {}
        for key, s in layer.items():
            for (kind, s2) in successors(g, s):
                s2 = closure(g, *s2)
                k2 = canon(g, *s2)
                if k2 not in new:
                    new[k2] = s2
                    if keep_examples: parent[k2] = (key, kind)
            if max_states and len(new) > max_states: break
        layer = new
        mx = max(maxlv(g, a) + maxlv(g, b) for (a, b, L) in layer.values())
        goal = [s for s in layer.values() if (s[0] & g.tmask) and (s[1] & g.tmask)]
        bad = [k for k, (a, b, L) in layer.items() if maxlv(g, a) + maxlv(g, b) > M]
        report("%s M=%d states=%d max(L1+L2)=%d (I) %s goal %s  [%.1fs]" % (
            g.name, M, len(layer), mx, "holds" if not bad else "VIOLATED in %d states" % len(bad),
            "REACHED (%d)" % len(goal) if goal else "no", time.time() - t0))
        hist.append((M, len(layer), mx, len(bad), len(goal)))
        if bad and keep_examples:
            k = bad[0]; chain = []
            while k in parent: chain.append((k, parent[k][1])); k = parent[k][0]
            report("   violation example (last state first):")
            for (kk, kind) in chain:
                a, b, L = kk
                report("     via %s: A1 lv %s, A2 lv %s, lineages %d" % (
                    kind, sorted(g.lv[x] for x in bits(a)), sorted(g.lv[x] for x in bits(b)), len(L)))
        if max_states and len(layer) > max_states:
            report("   (state cap reached; stopping)"); break
    return hist

# ---------------------------------------------------------------------------------------------- game builders
def from_toylab(G, P, hreps, target=None, name=""):
    from toylab import FrameGame
    fg = FrameGame(G, P, hreps)
    N = fg.ncos
    c0 = int(fg.coset[0])
    # reorder positions so that the raw frame is 0
    order = [c0] + [c for c in range(N) if c != c0]
    pos = {c: i for i, c in enumerate(order)}
    nb = [sorted(pos[c2] for c2 in fg.adj[c] if c2 != c) for c in order]
    lv = [fg.dist[c] for c in order]
    # twists: right multiplication P r -> P r t
    perms = {}
    ident = tuple(range(N))
    tws = [ident]; seen = {ident: 0}
    for t in range(G.n):
        p = tuple(pos[int(fg.coset[G.m(fg.reps[c], t)])] for c in order)
        if p not in seen: seen[p] = len(tws); tws.append(p)
    free = [i for i, p in enumerate(tws) if p[0] == 0]
    m = max(lv)
    if target is None:
        target = next(x for x in range(N) if lv[x] == m)
    orb = {tws[a][target] for a in free}
    g = Game(N, nb, lv, orb, tws, free, name=name, mult_table=(len(tws) <= 2000))
    g.tmask = sum(1 << x for x in orb)
    g.fg, g.order, g.pos, g.Gt = fg, order, pos, seen      # for algebraic checks: coset -> position index
    g.elem_pos = [pos[int(fg.coset[e])] for e in range(G.n)]
    g.elem_tw = [seen[tuple(pos[int(fg.coset[G.m(fg.reps[c], t)])] for c in order)] for t in range(G.n)]
    return g

if __name__ == "__main__":
    pass

# ---------------------------------------------------------------------------------------------- faster BFS (session 2)
def bfs2(g, D, report=print, max_states=None, last_nodedup=True, split=None):
    """layered BFS with raw dedup + canonical dedup up to depth D-1; at depth D (if last_nodedup) every successor of
    every depth-(D-1) state is checked for (I) and the goal without storing it.  Returns history rows."""
    st = closure(g, 1, 1, ())
    layer = {canon(g, *st): st}
    hist = []
    for M in range(1, D + 1):
        t0 = time.time()
        final = (M == D and last_nodedup)
        new = {}; seen_raw = set(); nsucc = 0; mx = -1; nbad = 0; ngoal = 0; ex = None
        items = sorted(layer.items())
        if split is not None and M == split[0] + 1:
            items = [it for i, it in enumerate(items) if i % split[2] == split[1]]
            report("   split: expanding part %d/%d of layer %d (%d states)" % (split[1], split[2], split[0], len(items)))
        for key, s in items:
            for (kind, s2) in successors(g, s):
                if s2 in seen_raw: continue
                seen_raw.add(s2)
                s2 = closure(g, *s2)
                nsucc += 1
                a, b = maxlv(g, s2[0]), maxlv(g, s2[1])
                if a + b > mx: mx = a + b
                if a + b > M:
                    nbad += 1
                    if ex is None: ex = (key, kind, s2)
                if (s2[0] & g.tmask) and (s2[1] & g.tmask): ngoal += 1
                if not final:
                    k2 = canon(g, *s2)
                    if k2 not in new: new[k2] = s2
            if len(seen_raw) > 4000000: seen_raw = set()
            if max_states and len(new) > max_states: break
        report("%s M=%d %s=%d max(L1+L2)=%d (I) %s goal %s  [%.1fs]" % (
            g.name, M, "successors(no dedup)" if final else "states", nsucc if final else len(new), mx,
            "holds" if not nbad else "VIOLATED (%d)" % nbad, "REACHED (%d)" % ngoal if ngoal else "no", time.time() - t0))
        if ex is not None:
            key, kind, s2 = ex
            report("   violation: from A1 lv %s A2 lv %s lins %s via %s -> A1 lv %s A2 lv %s lins %s" % (
                sorted(g.lv[x] for x in bits(key[0])), sorted(g.lv[x] for x in bits(key[1])), key[2], kind,
                sorted(g.lv[x] for x in bits(s2[0])), sorted(g.lv[x] for x in bits(s2[1])), s2[2]))
        hist.append((M, len(new), mx, nbad, ngoal))
        if final: break
        layer = new
        if max_states and len(layer) > max_states:
            report("   (state cap reached at M=%d; stopping)" % M); break
    return hist

def hamming_game(k, wreath=False, name=None):
    """product of k block-1 frame games: S_3^k (or S_3 wr S_k) on {0,1,2}^k, raw frame 0^k, moves change one coordinate.
    A sub-game of the block-k engine single-term game (frames that are direct sums in a fixed hyperbolic splitting)."""
    from toylab import perm_group, cyc
    n = 3 * k
    gens = []
    for i in range(k):
        gens.append(cyc(n, [3 * i, 3 * i + 1])); gens.append(cyc(n, [3 * i, 3 * i + 1, 3 * i + 2]))
    if wreath and k >= 2:
        sw = list(range(n))
        for j in range(3): sw[j], sw[3 + j] = 3 + j, j
        gens.append(tuple(sw))
        gens.append(tuple((x + 3) % n for x in range(n)))
    G = perm_group(gens, "H(%d,3)" % k)
    # P = stabiliser of the set {0, 3, 6, ...} pointwise-as-blocks (each block's base point 3i fixed)
    base = {3 * i for i in range(k)}
    P = [e for e, p in enumerate(G.elems) if {p[b] for b in base} == base]
    hreps = []
    for i in range(k):
        hreps.append(G.idx[cyc(n, [3 * i, 3 * i + 1])])
        if wreath: break
    return from_toylab(G, P, hreps, name=name or ("H%d3%s" % (k, "w" if wreath else "")))

# ---------------------------------------------------------------------------------------------- two final layers
def _locked_comps(g, lins):
    """vertex -> component id, and per component: max level of input-1 / input-2 vertices"""
    N = g.N
    par = {}
    def f(a):
        while par.get(a, a) != a:
            par[a] = par.get(par[a], par[a]); a = par[a]
        return a
    for (t, xs) in lins:
        for x in bits(xs):
            u, w = x, N + g.tw[t][x]
            ru, rw = f(u), f(w)
            if ru != rw: par[ru] = rw
    hi = {}
    verts = set(par)
    for (t, xs) in lins:
        for x in bits(xs): verts.add(x); verts.add(N + g.tw[t][x])
    comp = {}
    for v in verts:
        r = f(v); comp[v] = r
        h = hi.setdefault(r, [-1, -1])
        if v < N: h[0] = max(h[0], g.lv[v])
        else: h[1] = max(h[1], g.lv[v - N])
    return comp, hi

def goal_in_one(g, s, target_level):
    """exact: does some move from s make both inputs available at target_level (in the target orbit)?
    Moves that cannot matter are skipped: a frame of input l becomes available only as an end of the new generator
    or as a vertex of a component the new generator touches."""
    A1, A2, lins = s
    N = g.N
    need = [l for l, A in ((0, A1), (1, A2)) if not (A & g.tmask)]
    if not need: return True
    comp, hi = _locked_comps(g, lins)
    def touches(l, ends):
        # ends: list of (side, position) of the new generator
        for (side, p) in ends:
            if side == l and g.lv[p] >= target_level: return True
            v = p if side == 0 else N + p
            if v in comp and hi[comp[v]][l] >= target_level: return True
        return False
    for (kind, s2, ends) in successors(g, s, with_ends=True):
        if not all(touches(l, ends) for l in need): continue
        c = closure(g, *s2)
        if (c[0] & g.tmask) and (c[1] & g.tmask): return True
    return False

def bfs3(g, D, report=print, split=None):
    """dedup to depth D-2; depth D-1 every successor checked for (I) and the goal; depth D: exact goal-in-one test
    (goal at D = both inputs at level m; for D = 2m - 1 this is exactly a 2 in 2m - 1)."""
    st = closure(g, 1, 1, ())
    layer = {canon(g, *st): st}
    for M in range(1, D - 1):
        t0 = time.time(); new = {}
        items = sorted(layer.items())
        if split is not None and M == split[0] + 1:
            items = [it for i, it in enumerate(items) if i % split[2] == split[1]]
            report("   split: expanding part %d/%d of layer %d (%d states)" % (split[1], split[2], split[0], len(items)))
        mx = -1
        for key, s in items:
            for (kind, s2) in successors(g, s):
                s2 = closure(g, *s2); k2 = canon(g, *s2)
                if k2 not in new:
                    new[k2] = s2; mx = max(mx, maxlv(g, s2[0]) + maxlv(g, s2[1]))
        layer = new
        report("%s M=%d states=%d max(L1+L2)=%d (I) %s  [%.1fs]" % (g.name, M, len(layer), mx, "holds" if mx <= M else "VIOLATED", time.time() - t0))
    t0 = time.time(); n1 = 0; nbad = 0; ngoal1 = 0; ntest = 0; ngoal = 0; seen = set(); ex = None
    for key, s in layer.items():
        for (kind, s2) in successors(g, s):
            if s2 in seen: continue
            seen.add(s2)
            if len(seen) > 3000000: seen = set()
            c = closure(g, *s2); n1 += 1
            a, b = maxlv(g, c[0]), maxlv(g, c[1])
            if a + b > D - 1: nbad += 1
            if (c[0] & g.tmask) and (c[1] & g.tmask): ngoal1 += 1
            # admissible filter for a goal at depth D: each input not yet at the target needs a frame of level >= m-1
            Q = [-1, -1]
            for x in bits(c[0]): Q[0] = max(Q[0], g.lv[x])
            for y in bits(c[1]): Q[1] = max(Q[1], g.lv[y])
            for (t, xs) in c[2]:
                for x in bits(xs):
                    Q[0] = max(Q[0], g.lv[x]); Q[1] = max(Q[1], g.lv[g.tw[t][x]])
            if Q[0] < g.m - 1 or Q[1] < g.m - 1: continue
            ntest += 1
            if goal_in_one(g, c, g.m):
                ngoal += 1
                if ex is None: ex = (key, kind, c)
    report("%s M=%d raw successors=%d (I) %s goal %s | M=%d: %d states passed the filter, goal-in-one %s  [%.1fs]" % (
        g.name, D - 1, n1, "holds" if not nbad else "VIOLATED (%d)" % nbad, "REACHED (%d)" % ngoal1 if ngoal1 else "no",
        D, ntest, "REACHED (%d)" % ngoal if ngoal else "no", time.time() - t0))
    if ex: report("   example: %s" % (ex,))
