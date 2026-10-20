"""v-r-riders / indep.py -- independent re-check of r-riders (shares NO code with r-riders' toy.py / replay.py / lb.py).

Semantics used here (point evaluation, not formal sums): on f digit columns the Walsh basis vector whose column labels
have multiplicities m = (m_1..m_7) (m_0 is irrelevant) sees dir(z) as the scalar i^N_z(m) and shift(w) as (-1)^N_w(m),
N_z(m) = sum_y m_y <y,z>.  Values depend on m mod 4 only, so a register is a (4^7 x 6) complex array: its coefficient on
each of the 6 inputs at every point.  "Right at every m in Z4^7" == "right for every number of digit columns".

1. Replays the words of r-riders/words.txt (parsed letter by letter) plus A-plain and the B-riders word re-typed here,
   counts paid moves (= 'move' letters), checks every target register at all 16384 points and that every other
   register is zero.
2. Hub local gate rule, written from global_k.py (op 1: x != y, x is not the open centre, cur[x] == cur[y];
   op 2: only one centre open; op 3: centre equals its source's image), with per-register frame tracking, and
   a semantic one-frame check of every register after every letter (each input column is c * i^<m, g - s_l>).
3. Lower bound: my own search of flip sets D (pairs and all sizes) over all point pairs with equal shift values,
   rank by numpy then exact rank of the witness modulo p = 1000000009 (p = 1 mod 4); integer programme solved by
   brute force over n_z in 0..8 (no scipy).
usage: python3 -B indep.py"""
import sys, itertools, ast
pass  # numpy/scipy: install with pip (the original run used a local copy)
sys.dont_write_bytecode = True
import numpy as np

def ip(y, z): return bin(y & z).count('1') & 1
CHI = {z: np.array([ip(y, z) for y in range(8)]) for z in range(8)}
E = (1, 2, 4)
Ffr = CHI[1] + CHI[2] + CHI[4]                      # popcount, as integers (reduce mod 4 when used)
PSI = {a: Ffr - CHI[a] for a in E}
INP = ['x1', 'x2', 'x4', 'y1', 'y2', 'y4']
COL = {s: k for k, s in enumerate(INP)}
START = {'x%d' % a: CHI[a] for a in E}; START.update({'y%d' % a: 0 * CHI[1] for a in E})
PTS = np.array(list(itertools.product(range(4), repeat=7)), dtype=np.int64)   # m_1..m_7
P = len(PTS)
I4 = np.array([1, 1j, -1, -1j])
def ev(q):        # i^<m, q> at every point, q a frame on y = 0..7 (q(0) is irrelevant)
    return I4[(PTS @ (np.asarray(q)[1:] % 4)) % 4]
DIRV = {z: ev(CHI[z]) for z in range(1, 8)}
SHV = {w: ev(2 * CHI[w]) for w in range(1, 8)}

def term(inp, local, coef=1):
    """array (P,6) for input inp sitting at LOCAL frame `local` (physical = local - start)."""
    A = np.zeros((P, 6), dtype=complex); A[:, COL[inp]] = coef * ev(np.asarray(local) - START[inp]); return A

def targets(v):
    T = {}
    for a in E:
        x, y = 'x%d' % a, 'y%d' % a
        T['X%d' % a] = term(y, Ffr, -1) if v == 'A' else term(x, Ffr)
        T['Y%d' % a] = term(y, PSI[a]) + term(x, PSI[a])
    if v == 'B':
        T['Z'] = sum(term('y%d' % a, Ffr) + term('x%d' % a, Ffr) for a in E)
    return T

class Illegal(Exception): pass

class M:
    def __init__(s, hub):
        s.hub = hub; s.R = {}; s.fr = {}; s.paid = 0; s.centre = None; s.crec = None; s.start = dict(START)
        for a in E:
            s.R['X%d' % a] = term('x%d' % a, CHI[a]); s.fr['X%d' % a] = CHI[a] % 4
            s.R['Y%d' % a] = term('y%d' % a, 0 * CHI[1]); s.fr['Y%d' % a] = 0 * CHI[1]
    def g(s, r):
        if r not in s.R: s.R[r] = np.zeros((P, 6), dtype=complex); s.fr[r] = None
        return s.R[r]
    def oneframe(s, r, starts=None):
        """semantic check: every nonzero input column of r is c * i^<m, g - s_l> for ONE local frame g."""
        starts = starts or s.start; A = s.g(r); gs = set()
        for l in INP:
            col = A[:, COL[l]]
            if not np.any(np.abs(col) > 1e-9): continue
            c0 = col[0]
            if abs(c0) < 1e-9: raise Illegal('%s: column %s not a single term' % (r, l))
            q = []
            for y in range(1, 8):
                e = np.zeros(7, dtype=np.int64); e[y - 1] = 1
                k = int(np.flatnonzero((PTS == e).all(1))[0]); ratio = col[k] / c0
                q.append(int(np.argmin(np.abs(I4 - ratio))))
            q = np.array([0] + q)
            if not np.allclose(col, c0 * ev(q)): raise Illegal('%s: column %s holds several frames' % (r, l))
            gs.add(tuple((q + starts[l]) % 4))
        if len(gs) > 1: raise Illegal('%s holds terms at %d local frames' % (r, len(gs)))
        return next(iter(gs)) if gs else None
    def do(s, L):
        op = L[0]
        if op == 'move':
            _, r, z = L; s.R[r] = s.g(r) * DIRV[z][:, None]; s.paid += 1
            if s.fr.get(r) is not None: s.fr[r] = (s.fr[r] + CHI[z]) % 4
            if s.crec and s.crec[0] == r: s.crec[2].append(DIRV[z])
        elif op == 'shift':
            _, r, w = L; s.R[r] = s.g(r) * SHV[w][:, None]
            if s.fr.get(r) is not None: s.fr[r] = (s.fr[r] + 2 * CHI[w]) % 4
            if s.crec and s.crec[0] == r: s.crec[2].append(SHV[w])
        elif op == 'add':
            _, d, src, c = L
            if s.hub:
                if d == src: raise Illegal('add into itself')
                if d == s.centre: raise Illegal('add into the open centre (hub op 1: x != center)')
                fd, fs = s.fr.get(d), s.fr.get(src)
                if fd is not None and fs is not None and not np.array_equal(fd, fs):
                    raise Illegal('add %s += %s: local frames differ' % (d, src))
            s.g(d); s.R[d] = s.R[d] + c * s.g(src)
            if s.hub and s.fr.get(d) is None: s.fr[d] = s.fr.get(src)
        elif op == 'copy':
            _, C, src = L
            if s.hub and s.centre is not None: raise Illegal('second centre open')
            if np.any(np.abs(s.g(C)) > 1e-9): raise Illegal('centre not empty')
            s.R[C] = s.g(src).copy(); s.fr[C] = s.fr.get(src); s.centre = C; s.crec = (C, src, [])
        elif op == 'erase':
            _, C = L
            if s.hub:
                if s.crec is None or s.crec[0] != C: raise Illegal('erase of a non-centre')
                img = s.g(s.crec[1]).copy()
                for f in s.crec[2]: img = img * f[:, None]
                if not np.allclose(img, s.g(C)): raise Illegal('erase: centre is not its source image')
            s.R[C] = np.zeros((P, 6), dtype=complex); s.fr[C] = None
            if s.crec and s.crec[0] == C: s.crec = None; s.centre = None
        elif op == 'setlin':
            if s.hub: raise Illegal('setlin: overwrite is not a hub op')
            _, d, combo = L; s.R[d] = sum(c * s.g(r) for r, c in combo)
        else: raise SystemExit('unknown letter %r' % (L,))
        if s.hub:
            for r in s.R: s.oneframe(r)
    def check(s, T):
        bad = [r for r in T if not np.allclose(s.g(r), T[r])]
        bad += [r + '(dirty)' for r in s.R if r not in T and np.any(np.abs(s.R[r]) > 1e-9)]
        return bad
    def boundary(s):
        new = {}
        for a in E:
            new['x%d' % a] = (2 * CHI[a] - Ffr) % 4; new['y%d' % a] = (-PSI[a]) % 4
        out = {}
        for r in s.R:
            try: out[r] = s.oneframe(r, new)
            except Illegal as e: out[r] = 'ILLEGAL: ' + str(e)
        return out

def parse_words(path):
    W = {}; cur = None
    for line in open(path):
        if line.startswith('=='): cur = line.strip('= \n'); W[cur] = []
        elif line.strip(): W[cur].append(ast.literal_eval(line.split(None, 1)[1]))
    return W

PORTS = {a: [d for d in E if d != a] for a in E}
def word_A_plain():
    w = []
    for a in E:
        for z in PORTS[a]: w.append(('move', 'X%d' % a, z))
        for z in PORTS[a]: w.append(('move', 'Y%d' % a, z))
    for a in E: w += [('copy', 'C', 'X%d' % a), ('move', 'C', a), ('shift', 'C', a), ('add', 'Y%d' % a, 'C', 1), ('erase', 'C')]
    w.append(('BOUNDARY',))
    for a in E: w += [('copy', 'C', 'Y%d' % a), ('move', 'C', a), ('add', 'X%d' % a, 'C', -1), ('erase', 'C')]
    return w

def run(name, word, v):
    for hub in (True, False):
        m = M(hub); note = ''
        try:
            for L in word:
                if L[0] == 'BOUNDARY':
                    b = m.boundary(); bad = {r: g for r, g in b.items() if isinstance(g, str)}
                    note = ' | boundary: ' + ('all registers one-frame under new routes' if not bad else str(bad))
                    continue
                m.do(L)
            bad = m.check(targets(v))
            res = 'EXACT at all %d points, scratch clean' % P if not bad else 'WRONG: %s' % bad
        except Illegal as e: res = 'ILLEGAL: %s' % e
        print('%-46s %-4s %-9s paid %2d  %s%s' % (name, v, 'hub-local' if hub else 'wide', m.paid, res, note))

def rank_mod(A, p=1000000009):
    A = [list(r) for r in A]; rk = 0
    for c in range(len(A[0])):
        piv = next((i for i in range(rk, len(A)) if A[i][c] % p), None)
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]; iv = pow(A[rk][c], p - 2, p); A[rk] = [x * iv % p for x in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][c] % p:
                f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk

def lower_bound(v, Tdict=None):
    p = 1000000009; Ip = next(pow(g, (p - 1) // 4, p) for g in range(2, 50) if pow(g, (p - 1) // 2, p) == p - 1)
    assert Ip * Ip % p == p - 1
    T = Tdict or targets(v); rows = sorted(T); TM = np.stack([T[r] for r in rows], axis=1)   # (P, rows, 6)
    W4 = 4 ** np.arange(6, -1, -1)
    assert np.array_equal(PTS @ W4, np.arange(P))
    CHm = np.array([[ip(y, z) for y in range(1, 8)] for z in range(1, 8)])
    best = {}
    for d in itertools.product(range(4), repeat=7):
        nd = CHm @ np.array(d)
        if np.any(nd % 2): continue
        D = frozenset(z for z, x in zip(range(1, 8), nd % 4) if x == 2)
        if not D: continue
        k2 = ((PTS + np.array(d)) % 4) @ W4
        rk = np.linalg.matrix_rank(TM - TM[k2], tol=1e-7); j = int(np.argmax(rk))
        if rk[j] > best.get(D, (-1,))[0]: best[D] = (int(rk[j]), j, int(k2[j]))
    # exact re-check modulo p of every witness, with the shift / flip conditions re-derived from the points
    for D, (r, j, k) in best.items():
        m1, m2 = PTS[j], PTS[k]
        N1, N2 = CHm @ m1, CHm @ m2
        assert np.all((N1 - N2) % 2 == 0) and frozenset(np.flatnonzero((N1 - N2) % 4 == 2) + 1) == D
        def modp(x):
            re, im = int(round(x.real)), int(round(x.imag)); assert abs(x - complex(re, im)) < 1e-9
            return (re + im * Ip) % p
        Dm = [[modp(x) for x in row] for row in (TM[j] - TM[k])]
        assert rank_mod(Dm) == r, ('mod-p rank differs', sorted(D))
    Ds = list(best); req = np.array([best[D][0] for D in Ds]); A = np.array([[z in D for z in range(1, 8)] for D in Ds], dtype=np.int64)
    # brute-force integer programme: every rank is <= 6, so n_z in 0..6 loses nothing (lowering n_z > 6 to 6 keeps
    # every constraint); enumerate n_1 in an outer loop to keep memory small
    assert req.max() <= 6
    rest = np.array(list(itertools.product(range(7), repeat=6)), dtype=np.int16)
    best_tot, best_n = 10 ** 6, None
    for n1 in range(7):
        g = np.concatenate([np.full((len(rest), 1), n1, dtype=np.int16), rest], axis=1)
        ok = np.all(g @ A.T.astype(np.int16) >= req, axis=1)
        if ok.any():
            t = g.sum(1).astype(np.int64); t[~ok] = 10 ** 6; j = int(np.argmin(t))
            if t[j] < best_tot: best_tot, best_n = int(t[j]), g[j].tolist()
    pairs = sorted((sorted(D), best[D][0]) for D in Ds if len(D) == 2)
    return best_tot, best_n, pairs, {s: max(best[D][0] for D in Ds if len(D) == s) for s in sorted({len(D) for D in Ds})}

if __name__ == '__main__':
    W = parse_words(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../stage-sharing-riders/words.txt'))
    print('== 1+2. replays (my semantics; hub-local = global_k.py ops 1-3 + one-frame check every letter)')
    for name, w in W.items():
        v = name.split()[0]; run('words.txt: ' + name.split(' under')[0], w, v)
        if 'riders + cleanup' in name:
            noclean = [L for k, L in enumerate(w) if not (18 <= k <= 29)]
            run('  same word, cleanup letters 18-29 removed', noclean, 'B')
    run('A plain (re-typed here, with boundary)', word_A_plain(), 'A')
    # the B rider word, re-typed from v-s-wide's rider_check.py description (not from words.txt)
    wb = [('move', 'Y4', 1), ('move', 'Y1', 2), ('add', 'X1', 'Y4', 1), ('add', 'X2', 'Y1', 1),
          ('move', 'X1', 2), ('move', 'X1', 4), ('move', 'X2', 1), ('move', 'X2', 4), ('move', 'Y4', 2), ('move', 'Y1', 4),
          ('move', 'Y2', 1), ('move', 'Y2', 4), ('move', 'X4', 1), ('add', 'X4', 'Y2', 1), ('move', 'X4', 2)]
    wb += [('add', 'Z', 'X%d' % a, 1) for a in E]
    for a, s in ((1, 4), (2, 1), (4, 2)): wb += [('copy', 'C', 'Y%d' % s), ('move', 'C', s), ('add', 'X%d' % a, 'C', -1), ('erase', 'C')]
    for a in E: wb += [('copy', 'C', 'X%d' % a), ('move', 'C', a), ('shift', 'C', a), ('add', 'Y%d' % a, 'C', 1), ('erase', 'C')]
    run('B riders + cleanup (re-typed here)', wb, 'B')
    # skeptic probe: Y-side cleanup instead of X-side (riders swept into Y by the scatter, then removed there)
    wy = wb[:15] + [('add', 'Z', 'X%d' % a, 1) for a in E]
    for a in E: wy += [('copy', 'C', 'X%d' % a), ('move', 'C', a), ('shift', 'C', a), ('add', 'Y%d' % a, 'C', 1), ('erase', 'C')]
    run('probe: riders, scatter first, no X cleanup', wy, 'B')
    print('\n== 3. lower bound (own search, mod-p exact witnesses, brute-force IP)')
    for v in ('B', 'A'):
        b, n, pairs, per = lower_bound(v)
        print('[%s] bound on paid moves of EVERY word: %d  (n_z = %s); max rank per |D|: %s' % (v, b, n, per))
        print('    pair ranks: %s' % pairs)
