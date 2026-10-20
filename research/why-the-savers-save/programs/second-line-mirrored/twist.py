#!/usr/bin/env python3
"""cx-mirror / twist.py (own code, stdlib).  The twisted mirrors: every choice of second lines that keeps the v-code
isomorphic to the u-code.

v_b = s + j + A u_b with A in O(9, F_2) (A^T A = I for the dot product).  Then v_a . v_b = u_a . u_b for all a, b
(the V-centre can be a second copy of the P112 code), junk orthogonality holds automatically, and the remaining
condition u_b . v_b = 0 reads  q_A(u_b) := u_b . A u_b = 1  for all 112 ports.  N[a][b] = g_a . u_b = 1 + u_b . A u_a.
(1) The space Q of quadratic forms on F_2^9 that vanish on all 112 ports (45 coefficients: 9 linear, 36 cross).
    q_A - q_I has linear part diag(A) + 1 and cross part (A + A^T) off the diagonal; it must lie in Q.
(2) If Q = 0: A symmetric with unit diagonal, A^2 = I, so A = I + M with M symmetric, zero diagonal, M^2 = 0, and N is
    symmetric; a drop-free x chain is then an independent set of the graph {N = 1}.  We sample such A (products of
    Eichler maps x -> x + w (z.x) + z (w.x), w, z even and orthogonal, then keep only the A that pass the condition)
    and report, per A, the largest drop-free chain inside each total's support (exact maximum independent set by
    branch and bound on 56 vertices)."""
import gzip, json, os, random, itertools
from fractions import Fraction as F
P = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz'))
d = json.load(gzip.open(P)); v = d['v']; R = d['R']; n = 2 * v + R; U = d['ports']
dot = lambda a, b: bin(a & b).count('1') & 1
bit = lambda x, i: (x >> i) & 1
out = {}
# (1) Q
rows = []
pairs = [(i, j) for i in range(9) for j in range(i + 1, 9)]
for u in U:
    rows.append([bit(u, i) for i in range(9)] + [bit(u, i) & bit(u, j) for i, j in pairs])
def nullspace(rows, ncol):
    """F_2 nullspace of the matrix with the given rows"""
    piv = []; M = [r[:] for r in rows]; rr = 0; pc = []
    for c in range(ncol):
        p = next((i for i in range(rr, len(M)) if M[i][c]), None)
        if p is None: continue
        M[rr], M[p] = M[p], M[rr]
        for i in range(len(M)):
            if i != rr and M[i][c]: M[i] = [a ^ b for a, b in zip(M[i], M[rr])]
        pc.append(c); rr += 1
    free = [c for c in range(ncol) if c not in pc]
    basis = []
    for f in free:
        x = [0] * ncol; x[f] = 1
        for i, c in enumerate(pc):
            if M[i][f]: x[c] = 1
        basis.append(x)
    return basis
Q = nullspace(rows, 45)
print('dim Q (quadratic forms vanishing on all 112 ports):', len(Q))
out['dim_Q'] = len(Q)
# supports of the totals
cont = [({r: F(1)} if r < v else {}) for r in range(n)]
for g in d['A']:
    if g[0] == 'out':
        src = dict(cont[g[2]])
        for t, a, b in g[3]:
            for k, val in src.items():
                nv = cont[t].get(k, 0) + F(a, b) * val
                if nv: cont[t][k] = nv
                else: cont[t].pop(k, None)
    else:
        for s_, a, b in g[3]:
            for k, val in list(cont[s_].items()):
                nv = cont[g[2]].get(k, 0) + F(a, b) * val
                if nv: cont[g[2]][k] = nv
                else: cont[g[2]].pop(k, None)
Sk = [sorted(cont[s_]) for _, s_, _ in sorted(d['ret'])]

def apply(A, x):   # A as list of 9 column images (ints)
    y = 0
    for i in range(9):
        if bit(x, i): y ^= A[i]
    return y
I9 = [1 << i for i in range(9)]
def eichler(w, z):
    return [(1 << i) ^ (w if dot(z, 1 << i) else 0) ^ (z if dot(w, 1 << i) else 0) for i in range(9)]
def compose(A, B):  # A after B
    return [apply(A, B[i]) for i in range(9)]
def orth(A):
    return all(dot(A[i], A[j]) == (1 if i == j else 0) for i in range(9) for j in range(9))
def ok(A):
    return all(dot(u, apply(A, u)) == 1 for u in U)
def mis(nodes, adjN):
    best = []
    def bb(Rr, Pp):
        nonlocal best
        if len(Rr) + len(Pp) <= len(best): return
        if not Pp: best = list(Rr); return
        x = Pp[0]
        bb(Rr + [x], [y for y in Pp[1:] if y not in adjN[x]])
        bb(Rr, Pp[1:])
    bb([], sorted(nodes, key=lambda x: -len(adjN[x] & set(nodes)))); return best
def chains(A):
    N = [[0 if a == b else 1 ^ dot(U[b], apply(A, U[a])) for b in range(v)] for a in range(v)]
    sym = all(N[a][b] == N[b][a] for a in range(v) for b in range(v))
    adjN = {a: {b for b in range(v) if N[a][b] or N[b][a]} for a in range(v)}
    return sym, [len(mis(S, adjN)) for S in Sk]
even = [w for w in range(1, 512) if bin(w).count('1') % 2 == 0]
rnd = random.Random(7)
s0, c0 = chains(I9)
print('natural mirror A = I: N symmetric', s0, '| longest drop-free chain per support', c0)
out['mirror_chains'] = c0
found = []; tried = 0
while tried < 20000 and len(found) < 40:
    tried += 1
    A = I9
    for _ in range(rnd.randint(1, 4)):
        w = rnd.choice(even); z = rnd.choice([z for z in even if dot(w, z) == 0 and z != w])
        A = compose(eichler(w, z), A)
    if A != I9 and orth(A) and ok(A) and A not in found:
        found.append(A)
print('Eichler products tried', tried, '| distinct A != I passing u.Au = 1 on every port:', len(found))
res = []
for A in found:
    sym, c = chains(A)
    res.append(dict(A=A, symmetric=sym, chains=c))
best = max(res, key=lambda r: min(r['chains'])) if res else None
if res:
    print('over these A: longest drop-free chain in the worst support: max', min(best['chains']),
          '| per-support maxima over all A', [max(r['chains'][k] for r in res) for k in range(8)],
          '| all N symmetric', all(r['symmetric'] for r in res))
out['twisted'] = dict(tried=tried, found=len(found), best=best,
                      per_support_max=[max(r['chains'][k] for r in res) for k in range(8)] if res else None)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'twist.json'), 'w'), indent=1)
