"""p-pick sanity probe: for each isotropic h in F_2^9, the P112 ports in C_h = {S : u_S.h = 0}; the number of
P112 octads (8 mutually orthogonal ports) lying wholly inside C_h, and the largest orthogonal clique inside C_h."""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), "../break"))
from common import *
from collections import Counter
P = list(range(8, 120))
orth = {S: set(T for T in P if T != S and G[S][T] == 0) for S in P}
# octads inside P (8 mutually orthogonal ports)
octs = []
def rec(cur, cand):
    if len(cur) == 8:
        octs.append(tuple(cur)); return
    for T in sorted(cand):
        if cur and T < cur[-1]: continue
        rec(cur + [T], cand & orth[T])
for S in P:
    rec([S], set(T for T in orth[S] if T > S))
octs = set(tuple(sorted(o)) for o in octs)
print('octads in P112:', len(octs))
def omega(Vs):
    Vs = set(Vs); best = [0]
    def r(size, cand):
        if size + len(cand) <= best[0]: return
        if not cand:
            best[0] = max(best[0], size); return
        T = min(cand)
        r(size + 1, cand & orth[T]); r(size, cand - {T})
    r(0, Vs); return best[0]
res = Counter()
for h in range(1, 512):
    if dot9(h, h): continue
    C = [S for S in P if not dot9(UI[S], h)]
    Cs = set(C)
    inside = sum(1 for o in octs if set(o) <= Cs)
    res[(len(C), inside, omega(C))] += 1
for k, c in sorted(res.items()): print('|C_h| %d  octads inside %d  omega %d : %d isotropic h' % (k + (c,)))
