import sys, collections
sys.path.insert(0, '.')
from xbase import Field
Fd = Field(5)
# direction z (coordinate int) <-> field element: coord(x) = sum_j Tr(x f_j) 2^j
c2f = {Fd.coord(x): x for x in range(1, 32)}
Zline = open(sys.argv[1]).read().split('\n')
for ln in Zline:
    if ln.startswith('Z 0 255'):
        Z = [0] + [int(h, 16) for h in ln.split()[3:]]
ks = [int(x) for x in sys.argv[2].split(',')]
pts = set(Fd.exp[k % 31] for k in ks)
pairs = [(x, y) for x in range(1, 32) for y in range(x + 1, 32) if Z[x] >> (y - 1) & 1 and Z[y] >> (x - 1) & 1]
print('pairs', len(pairs))
# classify by the field ratio / sum of the two directions
cs = collections.Counter(); cr = collections.Counter()
for x, y in pairs:
    fx, fy = c2f[x], c2f[y]
    cs[Fd.log[fx ^ fy]] += 1
    r = Fd.log[fx] - Fd.log[fy]; r %= 31; r = min(r, (31 - r) % 31)
    cr[r] += 1
print('by log(x+y):', sorted(cs.items()))
print('by +-log(x/y):', sorted(cr.items()))
deg = collections.Counter()
for x, y in pairs: deg[x] += 1; deg[y] += 1
print('degrees by log of direction:', sorted((Fd.log[c2f[x]], deg[x]) for x in range(1, 32)))
print('target logs', ks)
