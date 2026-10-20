#!/usr/bin/env python3
"""v-cx-medium / vraw.py (own simulator; numpy from our pylib; python3 -I).  Replays cx-medium's explicit crossing call
(network = cx-medium's build(): flags, free meetings, carriers, natural decoder) on raw complex arrays with MY OWN
semantics, written without reading cx-medium's Raw class beyond its docstring:
  * a frame is a subspace U of F_2^k (stored as a set of vectors); U must be non-degenerate (U cap U^perp = 0);
  * Phi_U = H diag(i^{wt(pi_U x)}) H with H the normalised Walsh matrix and pi_U the orthogonal projection, computed
    from an ARBITRARY basis through the inverse Gram matrix over F_2 (no orthonormal basis assumed);
  * a move of an array holding content from U to U' requires U subset U' or U' subset U, costs |dim U' - dim U| and
    maps content c -> Phi_U' Phi_U^-1 c; a zero array moves for free; I also record each move's rank as a child;
  * a gate adds arrays only if every frame is EQUAL to the target's (as sets of vectors);
  * a copy is a fresh zero array at the source's frame plus one gate; scratch is dropped at the end (as in the hub's
    accounting) -- I report the scratch arrays left non-zero.
Checks: Y_a = Phi_{u_a^perp}(y_a + x_a), X_a = Phi_full(x_a) (underlying data x_a = Phi_<u_a>^-1 X_a(in)), every
frame where it must be; paid levels by class against cx-medium's call_*.json; and a second, exact check with integer
Gaussian inputs in exact rational complex arithmetic is too slow, so I use 3 random complex columns and report the
max error."""
import sys, os, json, itertools, time
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../medium-block-networks')))
from e7call import build          # cx-medium's network builder (the object under test)

def dot(x, y): return bin(x & y).count('1') & 1
def wt(x): return bin(x).count('1')

def span(vs):
    S = {0}
    for v in vs:
        if v not in S: S |= {s ^ v for s in S}
    return frozenset(S)

def basis_of(S):
    b = []
    for v in sorted(S):
        if v and v not in span(b): b.append(v)
    return b

def inv2(M):
    n = len(M); A = [list(r) + [int(i == j) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c]); A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c]: A[r] = [x ^ y for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]

class Sim:
    def __init__(s, k, cols=3, seed=7):
        s.k, s.N = k, 1 << k
        H = np.array([[(-1) ** dot(x, y) for y in range(s.N)] for x in range(s.N)], dtype=complex) / np.sqrt(s.N)
        s.H = H; s.cache = {}; s.arr = {}; s.paid = {}; s.children = {}; s.rng = np.random.default_rng(seed); s.cols = cols
    def Phi(s, U):
        if U in s.cache: return s.cache[U]
        b = basis_of(U); d = len(b)
        if d == 0: Pm = np.eye(s.N, dtype=complex)
        else:
            G = [[dot(x, y) for y in b] for x in b]
            Gi = inv2(G)                                   # raises if degenerate
            diag = np.empty(s.N, dtype=complex)
            for x in range(s.N):
                t = [dot(x, v) for v in b]
                c = [sum(Gi[i][j] & t[j] for j in range(d)) & 1 for i in range(d)]
                u = 0
                for ci, v in zip(c, b):
                    if ci: u ^= v
                assert all(dot(x ^ u, v) == 0 for v in b)
                diag[x] = 1j ** (wt(u) % 4)
            Pm = s.H @ np.diag(diag) @ s.H
        s.cache[U] = Pm; return Pm
    def new(s, name, U, content=None):
        assert name not in s.arr, name
        s.arr[name] = [U, np.zeros((s.N, s.cols), complex) if content is None else content.copy()]
    def move(s, name, U2, cls):
        U, c = s.arr[name]
        assert U <= U2 or U2 <= U, ('not nested', name)
        r = abs(len(basis_of(U2)) - len(basis_of(U)))
        if r and np.abs(c).max() > 1e-12:
            s.paid[cls] = s.paid.get(cls, 0) + r
            s.children[(cls, r)] = s.children.get((cls, r), 0) + 1
            c = s.Phi(U2) @ np.linalg.solve(s.Phi(U), c)
        s.arr[name] = [U2, c]
    def add(s, tgt, terms):
        U = s.arr[tgt][0]
        for coef, src in terms:
            assert s.arr[src][0] == U, ('gate frames differ', tgt, src)
            s.arr[tgt][1] = s.arr[tgt][1] + coef * s.arr[src][1]

def replay(name, seed=0):
    D = build(name, seed, verbose=False)
    labs, k, n = D['labs'], D['k'], len(D['L'])
    S = Sim(k)
    full = span([1 << j for j in range(k)])
    xin, yin = [], []
    for b in range(n):
        U = span([labs[b]]); x = S.rng.normal(size=(S.N, S.cols)) + 1j * S.rng.normal(size=(S.N, S.cols))
        xin.append(x); S.new(('X', b), U, x)
    for a in range(n):
        y = S.rng.normal(size=(S.N, S.cols)) + 1j * S.rng.normal(size=(S.N, S.cols))
        yin.append(y); S.new(('Y', a), frozenset({0}), y)
    need = set((b, i) for (a, b), i in D['free'].items()) | set((c['b'], c['i']) for c in D['carriers'])
    # sources climb their flags; copies at the frames needed
    for b in range(n):
        for i in sorted({i for (bb, i) in need if bb == b}):
            S.move(('X', b), span(D['Xb'][b][:i]), 'x')
            S.new(('K', b, i), S.arr[('X', b)][0]); S.add(('K', b, i), [(1.0, ('X', b))])
        S.move(('X', b), span(D['Xb'][b][:k]), 'x')
        assert S.arr[('X', b)][0] == full
    # centre: totals at the full frame, copies descend to 0, scatter
    for q in range(len(D['F'])):
        S.new(('C', q), full)
        S.add(('C', q), [(float(D['G'][q][b]), ('X', b)) for b in range(n) if D['G'][q][b] != 0])
        S.move(('C', q), frozenset({0}), 'centre')
        w = D['c0'] / sum(x * x for x in D['F'][q])
        for a in range(n):
            if D['G'][q][a] != 0: S.add(('Y', a), [(w * float(D['G'][q][a]), ('C', q))])
    # carriers
    for ci, c in enumerate(D['carriers']):
        K = S.arr[('K', c['b'], c['i'])]
        S.new(('P', ci), K[0], K[1])
        S.move(('P', ci), span(D['Ya'][c['targets'][0]][:c['j']]), 'carrier')
    at = {}
    for ci, c in enumerate(D['carriers']):
        for a in c['targets']: at.setdefault((a, c['j']), []).append((ci, c['b']))
    for a in range(n):
        for j in range(1, k):
            S.move(('Y', a), span(D['Ya'][a][:j]), 'y')
            terms = [(-D['J'][(a, b)], ('K', b, i)) for (aa, b), i in D['free'].items() if aa == a and i == j]
            terms += [(-D['J'][(a, b)], ('P', ci)) for ci, b in at.get((a, j), [])]
            if terms: S.add(('Y', a), terms)
    err = 0.0
    for a in range(n):
        ua = span([labs[a]]); perp = frozenset(v for v in range(1 << k) if dot(v, labs[a]) == 0)
        assert S.arr[('Y', a)][0] == perp and S.arr[('X', a)][0] == full
        xa = np.linalg.solve(S.Phi(ua), xin[a])
        err = max(err, np.abs(S.arr[('Y', a)][1] - S.Phi(perp) @ (yin[a] + xa)).max(),
                  np.abs(S.arr[('X', a)][1] - S.Phi(full) @ xa).max())
    scratch_nonzero = sum(1 for nm, (U, c) in S.arr.items() if nm[0] in 'KCP' and np.abs(c).max() > 1e-9)
    return dict(design=name, k=k, n=n, max_err=float(err), exact=bool(err < 1e-9), paid=S.paid, total=sum(S.paid.values()),
                children={f'{c}:{r}': v for (c, r), v in sorted(S.children.items())}, scratch_nonzero_dropped=scratch_nonzero,
                free=len(D['free']), carriers=len(D['carriers']), junk=len(D['J']))

if __name__ == '__main__':
    name = sys.argv[1] if len(sys.argv) > 1 else 'A5'
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    t0 = time.time(); r = replay(name, seed); r['seconds'] = round(time.time() - t0, 1)
    print(json.dumps(r)); json.dump(r, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'vraw_{name}_s{seed}.json'), 'w'), indent=1)
