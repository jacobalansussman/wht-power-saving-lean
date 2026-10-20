import sys
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-relays'))
from patterns import hamming_e8_hub
from collections import Counter
L, G = hamming_e8_hub()
def info(P, tag):
    sub = [[G[a][b] for b in P] for a in P]
    offd = Counter(abs(sub[i][j]) for i in range(9) for j in range(9) if i != j)
    # F2 independence of u-labels: rank over F2 of the integer roots mod 2? use L//? Use roots halved mod 2 not simple; use rank of G sub mod 2
    def rank2(rows):
        basis = []
        for x in rows:
            for b in basis: x = min(x, x ^ b)
            if x: basis.append(x)
        return len(basis)
    r2 = rank2([sum((sub[i][j] % 2) << j for j in range(9)) for i in range(9)])
    # number of live lines orthogonal to all / each puncture etc: degree profile
    print(tag, P, 'offdiag |G| counts', dict(offd), 'F2 rank of (G mod 2) on P', r2)
info([39, 52, 57, 71, 76, 78, 103, 109, 114], 'builder P (geo numbering)')
info([8, 17, 39, 51, 53, 58, 63, 72, 114], 'hub word P mapped to geo')
def prof(P):
    live = [i for i in range(120) if i not in P]
    c = Counter(sum(1 for p in P if G[i][p] != 0) for i in live)
    # second invariant: for each punctured line, number of other punctured lines non-orthogonal
    d = Counter(sum(1 for q in P if q != p and G[p][q] != 0) for p in P)
    return sorted(c.items()), sorted(d.items())
print(prof([39, 52, 57, 71, 76, 78, 103, 109, 114]))
print(prof([8, 17, 39, 51, 53, 58, 63, 72, 114]))
# is G mod 2 / root mod 2 independent? compute rank over F2 of the root vectors /? (labels u): use L entries mod 2 after halving where all even
import itertools
def u(i):
    x = L[i]
    if all(t % 2 == 0 for t in x): return tuple((t // 2) % 2 for t in x)
    return tuple(t % 2 for t in x)
for P in ([39, 52, 57, 71, 76, 78, 103, 109, 114], [8, 17, 39, 51, 53, 58, 63, 72, 114]):
    print([u(i) for i in P])
