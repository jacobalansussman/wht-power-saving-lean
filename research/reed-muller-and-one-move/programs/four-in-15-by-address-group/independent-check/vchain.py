#!/usr/bin/env python3
"""vchain.py (v-w8-rmexact, skeptic): independent re-implementation, written from the statements in PROOF.md only
(no code shared with w8-rmexact).  Tuple matrices over F_2, forms as dicts of coefficients, own set partitions,
own Konig vertex-cover bound (subset enumeration instead of max flow), own rectangle-partition recursion.

Conventions (from Theorem A): h acts on directions by x -> h x (matrix from row bitmasks: A[i][j] = bit j of row i),
on Sym by S -> h S h^T; moves [z z^T]; copy h of an array has target h h^T and receives one virtual move in direction
h z per physical move z.  Form Q(x) = sum_{i<=j} q_ij x_i x_j; <Q, S> = sum_{i<=j} q_ij S_ij.

Usage:
  vchain.py classes                       placement classes (O_4-conjugacy) for orders 3, 5, 15
  vchain.py chain W:M ORDER [rows...]     all four tests on every split, one group per fixed-point-free class
  vchain.py one ROWS W:M SPLIT            verbose verdict on one split (rows like 6,13,10,3; split like 3,3,3,3,3)
  vchain.py dirs ROWS W DIRS              verdict for an explicit list of physical directions
"""
import itertools, sys
from functools import lru_cache

m = 4
DIRS = list(range(1, 1 << m))
PAIRS = [(i, j) for i in range(m) for j in range(i, m)]

def vec(x): return [(x >> i) & 1 for i in range(m)]
def from_rows(rows): return tuple(tuple((rows[i] >> j) & 1 for j in range(m)) for i in range(m))
def to_rows(A): return [sum(A[i][j] << j for j in range(m)) for i in range(m)]
def mul(A, B): return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(m)) % 2 for j in range(m)) for i in range(m))
def tp(A): return tuple(tuple(A[j][i] for j in range(m)) for i in range(m))
def act(A, x):
    v = vec(x); return sum((sum(A[i][j] * v[j] for j in range(m)) % 2) << i for i in range(m))
ID = tuple(tuple(int(i == j) for j in range(m)) for i in range(m))

def det2(A):
    M_ = [list(r) for r in A]; n = m
    for c in range(n):
        p = next((r for r in range(c, n) if M_[r][c]), None)
        if p is None: return 0
        M_[c], M_[p] = M_[p], M_[c]
        for r in range(n):
            if r != c and M_[r][c]: M_[r] = [(a + b) % 2 for a, b in zip(M_[r], M_[c])]
    return 1

# forms: index = 10-bit integer over PAIRS
NF = 1 << len(PAIRS)
def qval(Q, x):
    v = vec(x); return sum(v[i] * v[j] for b, (i, j) in enumerate(PAIRS) if (Q >> b) & 1) % 2
def pairing(Q, S): return sum(S[i][j] for b, (i, j) in enumerate(PAIRS) if (Q >> b) & 1) % 2
SUPP = [frozenset(x for x in DIRS if qval(Q, x)) for Q in range(NF)]
assert len(set(SUPP)) == NF

def cyclic(A):
    G = [ID]; B = A
    while B != ID: G.append(B); B = mul(A, B)
    return G

def set_partitions(s):
    s = list(s)
    if not s: yield []; return
    a = s[0]
    for p in set_partitions(s[1:]):
        yield [[a]] + p
        for i in range(len(p)): yield p[:i] + [[a] + p[i]] + p[i + 1:]

class Game:
    def __init__(self, H, W):
        self.H = H; self.W = W; self.n = len(H)
        self.targets = [mul(h, tp(h)) for h in H]
        self.nodd = [W * sum(pairing(Q, S) for S in self.targets) for Q in range(NF)]
        orb = {}; self.orbits = []
        for x in DIRS:
            if x in orb: continue
            o = sorted({act(h, x) for h in H}); self.orbits.append(o)
            for y in o: orb[y] = len(self.orbits) - 1
        self.orb = orb
        self.free = all(len(o) == self.n for o in self.orbits)
        # containment lists for unions used
        self.cont = {}
        self.parts = {}
        for Q in range(1, NF):
            P = sorted(SUPP[Q])
            self.parts[Q] = list(set_partitions(P)) if len(P) <= 4 else [[[x] for x in P]]
    def contain(self, U):
        U = frozenset(U)
        if U not in self.cont: self.cont[U] = [G for G in range(1, NF) if U <= SUPP[G]]
        return self.cont[U]
    def counts_from_split(self, split):
        assert self.free
        return {x: split[self.orb[x]] for x in DIRS}
    def counts_from_dirs(self, dirs):
        C = {x: 0 for x in DIRS}
        for z in dirs:
            for h in self.H: C[act(h, z)] += 1
        return C
    def verdict(self, C, verbose=False):
        sl = [None] + [sum(C[x] for x in SUPP[Q]) - self.nodd[Q] for Q in range(1, NF)]
        for Q in range(1, NF):
            if sl[Q] < 0: return ('T1 deletion', Q, sorted(SUPP[Q]), sl[Q])
        msc = {}
        def ms(U):
            U = frozenset(U)
            if U not in msc: msc[U] = min(sl[G] for G in self.contain(U))
            return msc[U]
        # T2: cover rows with product partitions
        for Q in range(1, NF):
            for part in self.parts[Q]:
                tot = 0
                for A in part:
                    cA = sum(C[x] for x in A)
                    for B in part:
                        cB = sum(C[y] for y in B)
                        tot += min(cA, cB, ms(set(A) | set(B)))
                if tot < sl[Q]: return ('T2 cover', Q, sorted(SUPP[Q]), [len(A) for A in part], 'sum', tot, 'slack', sl[Q])
        # forced zeros
        Z = set()
        tight = [G for G in range(1, NF) if sl[G] == 0 and sum(C[x] for x in SUPP[G]) > 0]
        for G in tight:
            for x in SUPP[G]:
                for y in SUPP[G]: Z.add((x, y))
        # T3: Konig bound by subset enumeration
        for Q in range(1, NF):
            if sl[Q] == 0: continue
            P = [x for x in sorted(SUPP[Q]) if C[x] > 0]
            k = len(P)
            adj = [sum(1 << j for j, y in enumerate(P) if (x, y) not in Z) for x in P]
            cost = [C[x] for x in P]
            nb = [0] * (1 << k); cs = [0] * (1 << k)
            for s in range(1, 1 << k):
                low = (s & -s).bit_length() - 1
                nb[s] = nb[s & (s - 1)] | adj[low]; cs[s] = cs[s & (s - 1)] + cost[low]
            full = (1 << k) - 1
            best = min(cs[A] + cs[nb[full ^ A]] for A in range(1 << k))
            if best < sl[Q]: return ('T3 E-pattern', Q, sorted(SUPP[Q]), 'cover', best, 'slack', sl[Q])
        # T4: rectangle partitions of the non-forced-zero cells, supports with <= 4 live points
        for Q in range(1, NF):
            if sl[Q] == 0: continue
            P = [x for x in sorted(SUPP[Q]) if C[x] > 0]
            if len(P) > 4: continue
            cells = [(x, y) for x in P for y in P if (x, y) not in Z]
            if not cells: return ('T4 all cells forced zero', Q, sorted(SUPP[Q]), 'slack', sl[Q])
            cid = {c: i for i, c in enumerate(cells)}
            subs = [s for r in range(1, len(P) + 1) for s in itertools.combinations(P, r)]
            rect = {}
            for R in subs:
                for Cc in subs:
                    mk = 0
                    for x in R:
                        for y in Cc:
                            if (x, y) in cid: mk |= 1 << cid[(x, y)]
                    if mk == 0: continue
                    u = min(sum(C[x] for x in R), sum(C[y] for y in Cc), ms(set(R) | set(Cc)))
                    if mk not in rect or u < rect[mk]: rect[mk] = u
            rl = list(rect.items())
            @lru_cache(maxsize=None)
            def best(rem):
                if rem == 0: return 0
                low = rem & -rem
                return min((u + best(rem & ~mk) for mk, u in rl if (mk & low) and (mk & rem) == mk), default=10 ** 9)
            b = best((1 << len(cells)) - 1)
            if b < sl[Q]: return ('T4 rectangles', Q, sorted(SUPP[Q]), 'cover', b, 'slack', sl[Q])
        return ('pass', len(tight), 'tight')

def splits(M, k):
    for cut in itertools.combinations(range(M + k - 1), k - 1):
        b = (-1,) + cut + (M + k - 1,)
        yield tuple(b[i + 1] - b[i] - 1 for i in range(k))

def all_GL():
    return [A for A in (from_rows(r) for r in itertools.product(range(16), repeat=m)) if det2(A)]

def classes(orders):
    GL = all_GL(); print('|GL| =', len(GL))
    O4 = [g for g in GL if mul(g, tp(g)) == ID]; print('|O4| =', len(O4))
    inv = {g: next(h for h in GL if mul(g, h) == ID) for g in O4}
    out = {}
    for n in orders:
        subs = {}
        for g in GL:
            c = cyclic(g)
            if len(c) == n: subs.setdefault(frozenset(c), g)
        seen = set(); cl = []
        for s in sorted(subs, key=lambda s: sorted(s)):
            if s in seen: continue
            orbit = {frozenset(mul(mul(o, h), inv[o]) for h in s) for o in O4}
            seen |= orbit
            g = subs[s]; fpf = all(act(g, x) != x for x in DIRS)
            cl.append((len(orbit), to_rows(g), fpf, orbit))
        out[n] = (len(subs), cl)
        print('order %d: %d subgroups, %d classes; fixed-point-free class sizes %s' % (
            n, len(subs), len(cl), sorted([c[0] for c in cl if c[2]], reverse=True)))
    return out

if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'classes':
        classes([int(v) for v in (sys.argv[2] if len(sys.argv) > 2 else '3,5,15').split(',')])
    elif cmd == 'one':
        H = cyclic(from_rows([int(v) for v in sys.argv[2].split(',')]))
        W, M = map(int, sys.argv[3].split(':'))
        g = Game(H, W); sp = tuple(int(v) for v in sys.argv[4].split(','))
        print('|H| =', len(H), 'orbits', [o[0] for o in g.orbits], 'split', sp, '->', g.verdict(g.counts_from_split(sp)))
    elif cmd == 'dirs':
        H = cyclic(from_rows([int(v) for v in sys.argv[2].split(',')]))
        W = int(sys.argv[3]); dl = [int(v) for v in sys.argv[4].split(',')]
        g = Game(H, W)
        print('|H| =', len(H), 'W', W, 'dirs', dl, '->', g.verdict(g.counts_from_dirs(dl)))
    elif cmd == 'chain':
        W, M = map(int, sys.argv[2].split(':'))
        gens = [[int(v) for v in s.split(',')] for s in sys.argv[3:]]
        for rows in gens:
            H = cyclic(from_rows(rows)); g = Game(H, W)
            assert g.free
            surv = []; tot = 0; why = {}
            for sp in splits(M, len(g.orbits)):
                tot += 1
                r = g.verdict(g.counts_from_split(sp))
                why[r[0]] = why.get(r[0], 0) + 1
                if r[0] == 'pass': surv.append(sp)
            print('gen %s |H|=%d orbit reps %s  %d in %d: %d splits, %d survive %s; first-failing test counts %s' % (
                rows, len(H), [o[0] for o in g.orbits], W, M, tot, len(surv), surv if len(surv) <= 130 else surv[:20], why), flush=True)
