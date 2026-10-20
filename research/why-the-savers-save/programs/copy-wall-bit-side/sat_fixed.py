"""sat_fixed.py (d-bitspan) -- F_2 decoder (XOR model, any B fitting the bit pattern, any rank <= K) with the read
vector of total j confined to a FIXED flat: the live lines orthogonal to the listed lines (empty list = free).
Input: a match_*.json family from match_sat.py.  Gauge: none (flats differ per total).
usage: venv/bin/python -I sat_fixed.py out/match_s2_6x8_7x2.json"""
import sys, json, time; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from pysat.solvers import Cadical153
from pysat.formula import IDPool
from e8geo import *
d = json.load(open(sys.argv[1]))
P = {'s2': PUNCT_S2, 'word': PUNCT_WORD}[d['unit']]
live = [i for i in range(NL) if i not in P]; v = len(live)
fam = d['fam']; K = len(fam)
pool = IDPool(); cl = []
allowed = [[all(G[T, l] == 0 for l in f['lines']) for T in live] for f in fam]
X = [[pool.id(('x', j, t)) if allowed[j][t] else None for t in range(v)] for j in range(K)]
Y = [[pool.id(('y', j, s)) for s in range(v)] for j in range(K)]
def AND(a, b):
    z = pool.id(('and', a, b)); cl.extend([[-z, a], [-z, b], [z, -a, -b]]); return z
def XORc(lits, rhs):
    if not lits:
        if rhs: cl.append([])
        return
    acc = lits[0]
    for l in lits[1:]:
        z = pool.id(('xor', acc, l)); cl.extend([[-z, acc, l], [-z, -acc, -l], [z, -acc, l], [z, acc, -l]]); acc = z
    cl.append([acc] if rhs else [-acc])
for s in range(v):
    for t in range(v):
        if s == t or G[live[s], live[t]] != 0:
            XORc([AND(Y[j][s], X[j][t]) for j in range(K) if X[j][t] is not None], 1 if s == t else 0)
print('fixed flats', [f['lines'] for f in fam], 'lines allowed per total', [sum(a) for a in allowed], 'vars', pool.top, 'clauses', len(cl), flush=True)
sv = Cadical153(bootstrap_with=cl); t0 = time.time(); r = sv.solve()
print('RESULT', 'SAT' if r else 'UNSAT', '%.1f s' % (time.time() - t0), flush=True)
