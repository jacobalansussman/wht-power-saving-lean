"""probe_h2.py (g-break, H2) -- bit-type pattern (B[S,T] = 0 on non-orthogonal pairs, B[S,S] != 0): compare the
F_2 rank of the natural decoder (1 + G) mod 2 (when G is integral with G[S,S] even and off-diagonal in {0,+-1})
with the clique lower bound alpha = largest set of pairwise non-orthogonal lines, on root systems D4, D5, E6, E7, E8.
If F_2 rank > alpha somewhere small, that is where another payload field could beat F_2."""
import sys, itertools
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), "."))
from common import *

def lines_of(roots):
    out, seen = [], set()
    for r in roots:
        k = tuple(r); m = tuple(-x for x in r)
        if k in seen or m in seen: continue
        seen.add(k); out.append(list(r))
    return out

def Dn(n):
    rs = []
    for i, j in itertools.combinations(range(n), 2):
        for s in (1, -1):
            v = [0] * n; v[i] = 1; v[j] = s; rs.append(v)
    return lines_of(rs)

def sub_E8(cond):
    return [i for i in range(120) if cond(i)]

def f2rank(M):
    return f2dim([sum((M[i][j] & 1) << j for j in range(len(M))) for i in range(len(M))])

def alpha_nonorth(Gm):
    n = len(Gm)
    adj = [sum(1 << j for j in range(n) if j != i and Gm[i][j] != 0) for i in range(n)]
    best = [0]
    def rec(cur, cand):
        if cur + bin(cand).count('1') <= best[0]: return
        if not cand:
            best[0] = max(best[0], cur); return
        i = cand.bit_length() - 1
        rec(cur + 1, cand & adj[i])
        rec(cur, cand & ~(1 << i))
    rec(0, (1 << n) - 1)
    return best[0]

def report(name, Gm):
    n = len(Gm)
    B = [[(1 + Gm[i][j]) % 2 if i != j else 1 for j in range(n)] for i in range(n)]
    ok = all(B[i][j] == 0 for i in range(n) for j in range(n) if i != j and Gm[i][j] != 0)
    print('%-6s v %3d  natural (1+G) mod 2 fits: %s  F2 rank %d   alpha (pairwise non-orth) %d' %
          (name, n, ok, f2rank(B), alpha_nonorth(Gm)))

for n in (4, 5, 6):
    Ls = Dn(n)
    Gm = [[sum(a * b for a, b in zip(x, y)) for y in Ls] for x in Ls]
    report('D%d' % n, Gm)
# E7 = lines orthogonal to line 0 inside E8; E6 = lines orthogonal to an A2 (lines 0 and a line at 60 deg to it)
e7 = [i for i in range(120) if i != 0 and G[i][0] == 0]
report('E7', [[G[i][j] for j in e7] for i in e7])
j60 = next(j for j in range(120) if G[0][j] != 0 and j != 0)
e6 = [i for i in range(120) if i not in (0, j60) and G[i][0] == 0 and G[i][j60] == 0]
report('E6', [[G[i][j] for j in e6] for i in e6])
report('E8', [[G[i][j] for j in range(120)] for i in range(120)])
