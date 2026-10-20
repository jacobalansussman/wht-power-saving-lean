"""sat_light.py (d-bitspan) -- K = 10 F_2 decoders on the bit P111 pattern whose read span W contains `a` linearly
independent LIGHT vectors, each supported on the live lines of a root-orthogonal flat (class 6: (l,m)^perp, E6/D6;
class 7: l^perp, E7; class 5: rank-3 (l,m,n)^perp).  W is gauge invariant, so the clique gauge X|_C = [I_8; 0],
Y|_C = [I_8; 0] stays valid.  Light vector i = sum_j c_ij x_j; independence via D (10 x a) with c D = I_a.
Options: --caseB adds 'W does not contain W0' (z in W^perp, z not in W0^perp).
usage: venv/bin/python -I sat_light.py <s2|word> <classes e.g. 6x8> [--caseB] [--K 10]"""
import sys, time, itertools; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import numpy as np
from pysat.solvers import Cadical153
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
from e8geo import *
unit = sys.argv[1]
classes = []
for part in sys.argv[2].split(','):
    c_, n_ = part.split('x'); classes += [int(c_)] * int(n_)
classes.sort(reverse=True)
K = int(sys.argv[sys.argv.index('--K') + 1]) if '--K' in sys.argv else 10
P = {'s2': PUNCT_S2, 'word': PUNCT_WORD}[unit]
live = [i for i in range(NL) if i not in P]; v = len(live)
U = f2_labels()
C = cliques(live, 8)[0]; cpos = {live.index(T): i for i, T in enumerate(C)}
pool = IDPool(); cl = []
TRUE = pool.id('TRUE'); cl.append([TRUE])
X = [[(TRUE if cpos[t] == j else -TRUE) if t in cpos else pool.id(('x', j, t)) for t in range(v)] for j in range(K)]
Y = [[(TRUE if cpos[s] == j else -TRUE) if s in cpos else pool.id(('y', j, s)) for s in range(v)] for j in range(K)]
def AND(a, b):
    if a == TRUE: return b
    if b == TRUE: return a
    if a == -TRUE or b == -TRUE: return -TRUE
    z = pool.id(('and', a, b)); cl.extend([[-z, a], [-z, b], [z, -a, -b]]); return z
def XORv(lits):
    lits = [l for l in lits if l != -TRUE]; k = sum(1 for l in lits if l == TRUE); lits = [l for l in lits if l != TRUE]
    if not lits: return TRUE if k & 1 else -TRUE
    acc = lits[0]
    for l in lits[1:]:
        z = pool.id(('xor', acc, l)); cl.extend([[-z, acc, l], [-z, -acc, -l], [z, -acc, l], [z, acc, -l]]); acc = z
    return -acc if k & 1 else acc
def fix(lit, val):
    cl.append([lit] if val else [-lit])
for s in range(v):
    for t in range(v):
        if s == t or G[live[s], live[t]] != 0:
            fix(XORv([AND(Y[j][s], X[j][t]) for j in range(K)]), s == t)
if '--caseB' in sys.argv:
    z = [pool.id(('z', t)) for t in range(v)]
    for j in range(K): fix(XORv([AND(z[t], X[j][t]) for t in range(v)]), 0)
    cl.append([XORv([z[t] for t in range(v) if U[live[t]][k]]) for k in range(9)])
a = len(classes)
_FC = {}
def FLATS(r):
    if r not in _FC:
        F = flats(r)
        _FC[r] = [bits(m) for m in F[r]]
    return _FC[r]
c = [[pool.id(('c', i, j)) for j in range(K)] for i in range(a)]
LV = [[XORv([AND(c[i][j], X[j][t]) for j in range(K)]) for t in range(v)] for i in range(a)]
sel = []
for i, cc in enumerate(classes):
    if cc == 1:                                   # class 1: the light vector is a single live line e_T
        S1 = [pool.id(('one', i, t)) for t in range(v)]
        cl.extend(CardEnc.equals(lits=S1, bound=1, vpool=pool, encoding=EncType.seqcounter).clauses)
        for t in range(v):
            cl.append([-LV[i][t], S1[t]]); cl.append([LV[i][t], -S1[t]])
        sel.append([S1]); continue
    if cc in (2, 3):                              # class 2/3: support inside an (exactly enumerated) flat of rank <= cc
        FL = FLATS(cc)
        S1 = [pool.id(('fl', i, k)) for k in range(len(FL))]
        cl.extend(CardEnc.equals(lits=S1, bound=1, vpool=pool, encoding=EncType.seqcounter).clauses)
        cont = {}
        for k, f in enumerate(FL):
            for T in f: cont.setdefault(T, []).append(S1[k])
        for t in range(v):
            cl.append([-LV[i][t]] + cont.get(live[t], []))
        sel.append([S1]); continue
    nsel = {8: 0, 7: 1, 6: 2, 5: 3}[cc]
    S = [[pool.id(('sel', i, b, l)) for l in range(NL)] for b in range(nsel)]; sel.append(S)
    for b in range(nsel):
        cl.extend(CardEnc.equals(lits=S[b], bound=1, vpool=pool, encoding=EncType.seqcounter).clauses)
        for t in range(v):
            T = live[t]
            cl.append([-LV[i][t]] + [S[b][l] for l in range(NL) if G[T, l] == 0 and l != T])
    for b in range(1, nsel):
        for p in range(NL):
            for r in range(p, NL): cl.append([-S[b][p], -S[b - 1][r]])
    if nsel == 3:
        for p, r, s_ in itertools.combinations(range(NL), 3):
            if np.linalg.matrix_rank(np.array([L[p], L[r], L[s_]], float)) < 3: cl.append([-S[0][p], -S[1][r], -S[2][s_]])
for i in range(a - 1):                       # same class: sorted by first selector (permutation symmetry of light vectors)
    if classes[i] == classes[i + 1] and classes[i] not in (8, 2, 3):
        A_, B_ = sel[i][0], sel[i + 1][0]
        for p in range(len(A_)):
            for r in range(p): cl.append([-A_[p], -B_[r]])
D = [[pool.id(('D', j, k)) for k in range(a)] for j in range(K)]
for i in range(a):
    for k in range(a):
        fix(XORv([AND(c[i][j], D[j][k]) for j in range(K)]), i == k)
print('sat_light unit %s K %d light classes %s caseB %s vars %d clauses %d' % (unit, K, classes, '--caseB' in sys.argv, pool.top, len(cl)), flush=True)
sv = Cadical153(bootstrap_with=cl); t0 = time.time(); r = sv.solve()
print('RESULT', 'SAT' if r else 'UNSAT', '%.1f s' % (time.time() - t0), flush=True)
if r:
    m = set(l for l in sv.get_model() if l > 0)
    val = lambda l: (l == TRUE) or (l != -TRUE and l > 0 and l in m) or (l < 0 and l != -TRUE and -l not in m)
    for i in range(a):
        ch = [[l for l in range(len(S)) if S[l] in m] for S in sel[i]]
        supp = [live[t] for t in range(v) if val(LV[i][t])]
        print(' light %d class %d flat perp %s reads %d lines, Q-rank %d' % (i, classes[i], ch, len(supp), int(np.linalg.matrix_rank(L[supp].astype(float))) if supp else 0), flush=True)
