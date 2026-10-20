#!/usr/bin/env python3
"""copyfloor.py -- an exact small model of the copy wall (standard library only; run with python3 -I).

THE MODEL (decoder level).  A unit has v ports.  Port T has a frame label a_T, a vector of a frame space K^h
(K = Q for bit-type units, F_q for complex-type units), and a symmetric form N on K^h.  Source T starts at the
frame <a_T>; target S lives at frames inside a_S^perp.  Target S KNOWS source T (side information, cancellable
at no copy cost) iff T != S and N(a_S, a_T) = 0.

A DECODER over the payload field F_p is a list of terms (r_j, g_j), r_j, g_j in F_p^v.  Term j is one scratch
total r_j . x, delivered with weight g_j[S] into target S.  The delivered matrix B = sum_j g_j r_j^T must FIT:
B[S,S] != 0 and B[S,T] = 0 whenever T != S and S does not know T (N(a_S,a_T) != 0).

Costs of a term with read set R = supp r_j and served set D = supp g_j (P = span{a_T : T in R}):
  star  = dim P                                   (the total, formed at P, is dropped to frame 0: hub "beta_j >= k w_in")
  pair  = rank of the form N between P and span{a_S : S in D}
        = dim P - dim(P cap (span a_D)^perp)      (drop only to the common cut: #446 Theorem C / Gram level)
cst of a decoder = sum of term costs.  The STAR FLOOR is min over fitting decoders of the star cost; the PAIRING
FLOOR the same with pair.  pair <= star term by term, so PAIRING FLOOR <= STAR FLOOR <= v (one drop per port).

Commands (see README.md):
  python3 -I copyfloor.py e8                      calibrate on E8: bit P120 72, bit P111 63, complex P112 64
  python3 -I copyfloor.py floors <instance>       exact star and pairing floors by Dijkstra (small v only)
  python3 -I copyfloor.py eval <instance> <decoder.json>   check a candidate decoder: fits? star? pair?
  python3 -I copyfloor.py list                    built-in instances
"""
import heapq
import itertools
import json
import sys
from fractions import Fraction

# ------------------------------------------------------------------ linear algebra over Q or F_q
def rank(rows, field):
    """rank of a list of rows (lists of ints) over field ('Q' or a prime q)."""
    if field == 'Q':
        M = [[Fraction(x) for x in r] for r in rows]
    else:
        M = [[x % field for x in r] for r in rows]
    if not M:
        return 0
    nc = len(M[0]); rk = 0
    for c in range(nc):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        if field == 'Q':
            inv = 1 / M[rk][c]
        else:
            inv = pow(M[rk][c], field - 2, field)
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] * inv
                if field == 'Q':
                    M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
                else:
                    M[i] = [(a - f * b) % field for a, b in zip(M[i], M[rk])]
        rk += 1
        if rk == len(M):
            break
    return rk


def form(a, b, N, field):
    s = sum(a[i] * N[i][j] * b[j] for i in range(len(a)) for j in range(len(b)) if N[i][j])
    return s if field == 'Q' else s % field


class Unit:
    """v ports with frame labels over `field` (Q or prime q), form N, payload prime p."""
    def __init__(self, name, field, labels, N=None, p=2, note=''):
        self.name, self.field, self.labels, self.p, self.note = name, field, [list(a) for a in labels], p, note
        self.h = len(labels[0])
        self.N = N if N is not None else [[int(i == j) for j in range(self.h)] for i in range(self.h)]
        self.v = len(labels)
        self.G = [[form(a, b, self.N, field) for b in self.labels] for a in self.labels]
        for T in range(self.v):
            if self.G[T][T] == 0:
                raise ValueError('label %d is isotropic: a target could read its own source at no cost' % T)
        # knows[S][T]: S may cancel x_T (side information)
        self.knows = [[S != T and self.G[S][T] == 0 for T in range(self.v)] for S in range(self.v)]
        self._star, self._pair = {}, {}

    def star(self, R):
        key = frozenset(R)
        if key not in self._star:
            self._star[key] = rank([self.labels[T] for T in R], self.field)
        return self._star[key]

    def pair(self, R, D):
        key = (frozenset(R), frozenset(D))
        if key not in self._pair:
            self._pair[key] = rank([[self.G[S][T] for T in sorted(R)] for S in sorted(D)], self.field) if R and D else 0
        return self._pair[key]

    # -------------------------------------------------------------- decoders
    def delivered(self, terms):
        p = self.p
        B = [[0] * self.v for _ in range(self.v)]
        for r, g in terms:
            for S in range(self.v):
                if g[S] % p:
                    for T in range(self.v):
                        if r[T] % p:
                            B[S][T] = (B[S][T] + g[S] * r[T]) % p
        return B

    def fits(self, B):
        bad = []
        for S in range(self.v):
            if B[S][S] % self.p == 0:
                bad.append(('diagonal zero', S))
            for T in range(self.v):
                if T != S and not self.knows[S][T] and B[S][T] % self.p:
                    bad.append(('interference S does not know', S, T))
        return bad

    def cost(self, terms):
        st = pr = 0
        rows = []
        for r, g in terms:
            R = [T for T in range(self.v) if r[T] % self.p]
            D = [S for S in range(self.v) if g[S] % self.p]
            s, q = self.star(R), self.pair(R, D)
            st += s; pr += q
            rows.append({'read': R, 'served': D, 'star': s, 'pair': q})
        return st, pr, rows


# ------------------------------------------------------------------ exact floors by Dijkstra (small v)
def exact_floor(U, kind='star', limit_states=5_000_000, verbose=False):
    """min over fitting decoders of the star (or pair) cost, exactly.  State = values of B on the constrained
    entries (diagonal and forbidden pairs); a move adds one term.  Terms that touch no constrained entry are
    useless and skipped.  For each effect on the constrained entries only the cheapest term is kept."""
    p, v = U.p, U.v
    E = [(S, S) for S in range(v)] + [(S, T) for S in range(v) for T in range(v) if S != T and not U.knows[S][T]]
    idx = {e: i for i, e in enumerate(E)}
    nE = len(E)
    # nonzero vectors up to scaling for r, all nonzero for g (scaling absorbed)
    def vecs(normalise):
        out = []
        for t in itertools.product(range(p), repeat=v):
            if any(t):
                if normalise and t[next(i for i in range(v) if t[i])] != 1:
                    continue
                out.append(t)
        return out
    R_list, G_list = vecs(True), vecs(False)
    best = {}
    for r in R_list:
        Rs = [T for T in range(v) if r[T]]
        for g in G_list:
            Ds = [S for S in range(v) if g[S]]
            eff = [0] * nE
            nz = False
            for S in Ds:
                for T in Rs:
                    i = idx.get((S, T))
                    if i is not None:
                        eff[i] = g[S] * r[T] % p
                        nz = nz or eff[i] != 0
            if not nz:
                continue
            c = U.star(Rs) if kind == 'star' else U.pair(Rs, Ds)
            key = tuple(eff)
            if key not in best or c < best[key][0]:
                best[key] = (c, (r, g))
    moves = [(c, eff, tg) for eff, (c, tg) in best.items()]
    moves.sort(key=lambda m: m[0])
    def is_goal(st):
        return all(st[i] % p for i in range(v)) and all(st[i] % p == 0 for i in range(v, nE))
    start = tuple([0] * nE)
    dist = {start: 0}; prev = {}
    pq = [(0, start)]
    n = 0
    while pq:
        d, st = heapq.heappop(pq)
        if d > dist.get(st, 1 << 30):
            continue
        if is_goal(st):
            path = []
            while st in prev:
                st, tg = prev[st]
                path.append(tg)
            return d, path[::-1], n
        n += 1
        if n > limit_states:
            return None, None, n
        for c, eff, tg in moves:
            ns = tuple((a + b) % p for a, b in zip(st, eff))
            nd = d + c
            if nd < dist.get(ns, 1 << 30):
                dist[ns] = nd; prev[ns] = (st, tg)
                heapq.heappush(pq, (nd, ns))
    return None, None, n


# ------------------------------------------------------------------ fast exact STAR floor over F_2 (branch and bound)
def star_floor_f2(U, verbose=False):
    """Exact star floor for payload F_2.  The star cost depends only on the read vectors r_j, and a family of read
    vectors admits served vectors making B fit iff for every target S the span W = span{r_j} contains b with
    b_S = 1 and b_T = 0 for every T that S does not know.  So the star floor is the least total weight of an
    independent family of read vectors whose span passes that test (weight w(r) = dim span of the labels read).
    Depth-first search over read vectors in (weight, mask) order, with branch and bound.  Returns (floor, family)."""
    assert U.p == 2
    v = U.v
    forb = [sum(1 << T for T in range(v) if T != S and not U.knows[S][T]) for S in range(v)]
    w = {}
    for r in range(1, 1 << v):
        w[r] = U.star([T for T in range(v) if r >> T & 1])
    order = sorted(w, key=lambda r: (w[r], r))
    def in_span(x, basis):
        for b in basis:
            x = min(x, x ^ b)
        return x == 0
    def valid(fam):
        for S in range(v):
            M = forb[S] | (1 << S)
            proj = []
            for r in fam:
                x = r & M
                for b in proj:
                    x = min(x, x ^ b)
                if x:
                    proj.append(x); proj.sort(reverse=True)
            if not in_span(1 << S, proj):
                return False
        return True
    best = [v + 1, None]
    # upper bound: the trivial decoder (one term per port) has star cost = sum of 1 = v
    best = [sum(w[1 << T] for T in range(v)), [1 << T for T in range(v)]]
    nodes = [0]
    def dfs(start, fam, basis, cost):
        nodes[0] += 1
        if fam and valid(fam):
            if cost < best[0]:
                best[0], best[1] = cost, list(fam)
            return                      # adding more vectors only costs more
        for k in range(start, len(order)):
            r = order[k]
            if cost + w[r] >= best[0]:
                break                   # order is by weight: nothing later is cheaper
            if in_span(r, basis):
                continue
            nb = sorted(basis + [r], reverse=True)
            # keep an echelon-ish basis for in_span: reduce
            red = []
            for b in nb:
                x = b
                for c in red:
                    x = min(x, x ^ c)
                if x:
                    red.append(x); red.sort(reverse=True)
            fam.append(r)
            dfs(k + 1, fam, red, cost + w[r])
            fam.pop()
    dfs(0, [], [], 0)
    return best[0], best[1], nodes[0]


# ------------------------------------------------------------------ E8 (Hamming construction, as in the hub's e8units.py)
def e8_lines():
    H = [0b00001111, 0b00110011, 0b01010101, 0b11111111]   # generators of the extended Hamming [8,4,4] code
    words = set()
    for c in itertools.product((0, 1), repeat=4):
        w = 0
        for bit, g in zip(c, H):
            if bit:
                w ^= g
        words.add(w)
    supports = sorted(tuple(i for i in range(8) if w >> i & 1) for w in words if bin(w).count('1') == 4)
    assert len(supports) == 14
    L = []
    for a in range(8):
        L.append([2 if i == a else 0 for i in range(8)])
    for sup in supports:
        for signs in itertools.product((1, -1), repeat=3):
            vec = [0] * 8
            vec[sup[0]] = 1
            for s, i in zip(signs, sup[1:]):
                vec[i] = s
            L.append(vec)
    assert len(L) == 120
    G = [[sum(x * y for x, y in zip(a, b)) // 2 for b in L] for a in L]
    return L, G


def f2_labels(B):
    """U (rows as ints) with U U^T = B over F_2, B symmetric non-alternating; columns = rank(B)."""
    import random
    n = len(B)
    rows = [sum((B[i][j] & 1) << j for j in range(n)) for i in range(n)]
    # row basis
    basis, I = {}, []
    for i, x in enumerate(rows):
        y = x
        while y:
            hb = y.bit_length() - 1
            if hb in basis:
                y ^= basis[hb]
            else:
                basis[hb] = y; I.append(i); break
    r = len(I)
    M = [[B[i][j] & 1 for j in I] for i in I]          # r x r invertible symmetric
    def dot(x, y):
        return sum(x[i] * M[i][j] * y[j] for i in range(r) for j in range(r)) & 1
    rng = random.Random(1)
    for _ in range(2000):                              # randomised Gram-Schmidt for an M-orthonormal basis X
        X = []
        ok = True
        space = [[int(i == j) for j in range(r)] for i in range(r)]
        for k in range(r):
            # vectors in the M-orthogonal complement of X: search random combinations of a complement basis
            cands = []
            for _t in range(400):
                c = [0] * r
                for b in space:
                    if rng.random() < 0.5:
                        c = [(x + y) & 1 for x, y in zip(c, b)]
                if any(c) and dot(c, c) == 1 and all(dot(c, x) == 0 for x in X):
                    cands.append(c); break
            if not cands:
                ok = False; break
            X.append(cands[0])
            # recompute complement basis: solve dot(y, x) = 0 for all x in X
            cons = [[sum(M[i][j] * x[j] for j in range(r)) & 1 for i in range(r)] for x in X]
            space = nullspace_f2(cons, r)
        if ok and len(X) == r:
            break
    else:
        raise RuntimeError('no orthonormal basis found')
    # U = B[:, I] X   (X columns)
    U = []
    for i in range(n):
        row = [B[i][j] & 1 for j in I]
        U.append([sum(row[a] * X[c][a] for a in range(r)) & 1 for c in range(r)])
    for i in range(n):
        for j in range(n):
            assert sum(U[i][c] * U[j][c] for c in range(r)) & 1 == B[i][j] & 1
    return U


def nullspace_f2(cons, r):
    """basis of {y in F_2^r : c.y = 0 for c in cons}."""
    M = [list(c) for c in cons]
    piv_cols, rk = [], 0
    for c in range(r):
        p = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                M[i] = [(a + b) & 1 for a, b in zip(M[i], M[rk])]
        piv_cols.append(c); rk += 1
    free = [c for c in range(r) if c not in piv_cols]
    out = []
    for f in free:
        y = [0] * r; y[f] = 1
        for i, c in enumerate(piv_cols):
            y[c] = M[i][f]
        out.append(y)
    return out


def e8_calibration():
    L, G = e8_lines()
    B0 = [[(1 + G[i][j]) % 2 for j in range(120)] for i in range(120)]
    U = f2_labels(B0)
    out = {}
    # ---- bit P120: frames over Q (the 120 lines), decoder = the 9 columns of U (B0 = U U^T)
    def bit_unit(live):
        return Unit('e8-bit', 'Q', [L[i] for i in live], p=2)
    live = list(range(120))
    Ub = bit_unit(live)
    terms = [([U[i][c] for i in live], [U[i][c] for i in live]) for c in range(9)]
    B = Ub.delivered(terms)
    assert not Ub.fits(B), Ub.fits(B)[:3]
    st, pr, rows = Ub.cost(terms)
    out['bit P120, 9 label totals (columns of U)'] = (st, pr)
    # every nonzero functional q: reader span
    spans = []
    for q in range(1, 512):
        R = [i for i in live if sum(U[i][c] for c in range(9) if q >> c & 1) & 1]
        spans.append(Ub.star(R))
    out['bit P120: min reader span over all 511 functionals'] = min(spans)
    # ---- bit P111: hub's shipped puncture; greedy min-weight basis of functionals (matroid), with dual served sets
    # The hub's shipped puncture [39,52,57,71,76,78,103,109,114] is in the hub's own line numbering.  What matters
    # for the copy floor is that the nine punctured lines have F_2-independent labels u_y (then the functional
    # q = u_y reads only y and y-perp, which loses a dimension when y is removed).  We take the first such set
    # in our numbering, one line per 4-frame of mutually orthogonal lines is not needed for this count.
    P, basis = [], {}
    for i in range(8, 120):
        x = sum(U[i][c] << c for c in range(9))
        y = x
        while y:
            hb = y.bit_length() - 1
            if hb in basis:
                y ^= basis[hb]
            else:
                basis[hb] = y; P.append(i); break
        if len(P) == 9:
            break
    out['bit P111: punctured lines (our numbering)'] = P
    live = [i for i in range(120) if i not in P]
    Ub = bit_unit(live)
    U9 = [U[i] for i in live]
    def reader(q):
        return [k for k in range(len(live)) if sum(U9[k][c] for c in range(9) if q >> c & 1) & 1]
    w = {q: Ub.star(reader(q)) for q in range(1, 512)}
    basis, chosen = [], []
    for q in sorted(range(1, 512), key=lambda q: (w[q], q)):
        x = q
        for b in sorted(basis, reverse=True):
            x = min(x, x ^ b)
        if x:
            basis.append(x); chosen.append(q)
        if len(chosen) == 9:
            break
    # decoder: B = U9 U9^T = sum_j (U9 q_j)(U9 q*_j)^T with Q = [q_j], Q^{-T} = [q*_j]
    Qm = [[chosen[j] >> c & 1 for j in range(9)] for c in range(9)]    # 9x9, columns q_j
    Qinv = inv_f2(Qm)
    QinvT = [[Qinv[j][c] for j in range(9)] for c in range(9)]
    duals = [sum(QinvT[c][j] << c for c in range(9)) for j in range(9)]
    terms = []
    for q, qs in zip(chosen, duals):
        r = [sum(U9[k][c] for c in range(9) if q >> c & 1) & 1 for k in range(len(live))]
        g = [sum(U9[k][c] for c in range(9) if qs >> c & 1) & 1 for k in range(len(live))]
        terms.append((r, g))
    B = Ub.delivered(terms)
    assert not Ub.fits(B)
    st, pr, rows = Ub.cost(terms)
    out['bit P111 (independent puncture), greedy star basis'] = (st, pr)
    out['bit P111: star weights of the chosen totals'] = [w[q] for q in chosen]
    # ---- complex P112 (hub cxp.py formula): totals A_a = live lines non-orthogonal to f_a = 2e_a; star = F2 span of
    # the 9-bit labels u_A; ThmC = F2 rank of B0[A, A].  Reproduced as the hub computes it.
    live = list(range(8, 120))
    st = th = 0
    for a in range(8):
        A = [S for S in live if G[S][a] != 0]
        st += rank([U[S] for S in A], 2)
        th += rank([[B0[i][j] for j in A] for i in A], 2)
    out['complex P112 frame decoder (hub formula): star, ThmC'] = (st, th)
    # ---- E7 (63 lines of E8 orthogonal to line 0), bit-type: smallest root system with x < 1 (hub ade_screen)
    e7 = [i for i in range(120) if G[i][0] == 0 and i != 0]
    B7 = [[B0[i][j] for j in e7] for i in e7]
    U7 = f2_labels(B7)
    c7 = len(U7[0])
    U7u = Unit('e7-bit', 'Q', [L[i] for i in e7], p=2)
    w7 = {}
    for q in range(1, 1 << c7):
        R = [k for k in range(63) if sum(U7[k][c] for c in range(c7) if q >> c & 1) & 1]
        w7[q] = U7u.star(R)
    basis, chosen = [], []
    for q in sorted(w7, key=lambda q: (w7[q], q)):
        x = q
        for b in sorted(basis, reverse=True):
            x = min(x, x ^ b)
        if x:
            basis.append(x); chosen.append(q)
    out['E7 bit: v, F2 rank of the pattern B0, label-decoder star floor (greedy), weights'] = (
        63, c7, sum(w7[q] for q in chosen), sorted(w7[q] for q in chosen))
    return out


def inv_f2(M):
    n = len(M)
    A = [list(M[i]) + [int(i == j) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if A[i][c])
        A[c], A[p] = A[p], A[c]
        for i in range(n):
            if i != c and A[i][c]:
                A[i] = [(x + y) & 1 for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


# ------------------------------------------------------------------ built-in small instances
def odd_vectors(h):
    return [[x >> i & 1 for i in range(h)] for x in range(1, 1 << h) if bin(x).count('1') % 2]


def instances():
    I = {}
    I['f2h3-odd'] = dict(field=2, labels=odd_vectors(3), note='complex-type: all 4 non-isotropic vectors of F_2^3')
    I['f2h4-odd'] = dict(field=2, labels=odd_vectors(4), note='complex-type: all 8 non-isotropic vectors of F_2^4')
    # A2 and A3 = D3 root lines (bit-type, frames over Q)
    I['A2'] = dict(field='Q', labels=[[1, -1, 0], [0, 1, -1], [1, 0, -1]], note='bit-type: 3 lines, pairwise non-orthogonal')
    I['A3'] = dict(field='Q', labels=[[1, -1, 0, 0], [1, 0, -1, 0], [1, 0, 0, -1], [0, 1, -1, 0], [0, 1, 0, -1],
                                      [0, 0, 1, -1]], note='bit-type: the 6 root lines of A3 = D3')
    I['B2'] = dict(field='Q', labels=[[1, 0], [0, 1], [1, 1], [1, -1]], note='bit-type: 4 lines of B2 (two orthogonal pairs)')
    I['cube13'] = dict(field='Q', labels=[[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1],
                                          [0, 1, 1], [0, 1, -1], [1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1]],
                       note='bit-type: the 13 lines of the cube (axes, face and body diagonals)')
    I['D4'] = dict(field='Q', labels=[list(v) for v in sorted({tuple(x) for x in _d4()})], note='bit-type: 12 D4 root lines')
    return I


def _d4():
    out = []
    for i, j in itertools.combinations(range(4), 2):
        for s in (1, -1):
            v = [0] * 4; v[i] = 1; v[j] = s; out.append(v)
    return out


def make(name, p=2, dup=1):
    spec = instances()[name]
    labels = [a for a in spec['labels'] for _ in range(dup)]
    return Unit(name + ('' if dup == 1 else '-x%d' % dup), spec['field'], labels, p=p, note=spec['note'])


# ------------------------------------------------------------------ CLI
def main(argv):
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__); return 0
    cmd = argv[0]
    if cmd == 'list':
        for k, s in instances().items():
            print('%-10s field %-3s v %3d  %s' % (k, s['field'], len(s['labels']), s['note']))
        return 0
    if cmd == 'e8':
        for k, val in e8_calibration().items():
            print('%-60s %s' % (k, val))
        return 0
    if cmd == 'floors':
        name = argv[1]
        p = int(argv[argv.index('--p') + 1]) if '--p' in argv else 2
        dup = int(argv[argv.index('--dup') + 1]) if '--dup' in argv else 1
        U = make(name, p, dup)
        print('instance %s: v %d, h %d, frame field %s, payload F_%d; knows-pairs %d' %
              (U.name, U.v, U.h, U.field, U.p, sum(map(sum, U.knows))))
        if p == 2 and '--dijkstra' not in argv:
            val, fam, n = star_floor_f2(U)
            print('  star floor = %d   (branch and bound, nodes %d); read sets of an optimal decoder:' % (val, n))
            for r in fam:
                R = [T for T in range(U.v) if r >> T & 1]
                print('     read %s  star %d' % (R, U.star(R)))
            kinds = ('pair',) if U.v <= 5 or '--pair' in argv else ()
        else:
            kinds = ('star', 'pair')
        for kind in kinds:
            val, path, n = exact_floor(U, kind)
            print('  %s floor = %s   (states expanded %d)' % (kind, val, n))
            if path:
                st, pr, rows = U.cost(path)
                for r in rows:
                    print('     term read %s served %s  star %d pair %d' % (r['read'], r['served'], r['star'], r['pair']))
        print('  one drop per port = %d' % U.v)
        return 0
    if cmd == 'eval':
        name, path = argv[1], argv[2]
        d = json.load(open(path))
        U = make(name, d.get('p', 2), d.get('dup', 1))
        terms = [(t['r'], t['g']) for t in d['terms']]
        B = U.delivered(terms)
        bad = U.fits(B)
        st, pr, rows = U.cost(terms)
        print('fits: %s%s' % (not bad, '' if not bad else '  first problems: %s' % bad[:5]))
        print('star cost %d, pairing cost %d, ports %d' % (st, pr, U.v))
        return 0 if not bad else 1
    print(__doc__); return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
