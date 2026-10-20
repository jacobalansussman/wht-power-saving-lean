# thmE.py: check every number in Theorem E with own GF(16) (x^4 + x + 1, polynomial basis, bit i = coeff of alpha^i)
import sys, itertools
sys.argv = ['x']
exec(open('vchain.py').read().split("if __name__ == '__main__':")[0])
def gmul(a, b):
    r = 0
    for i in range(4):
        if (b >> i) & 1: r ^= a << i
    for d in (6, 5, 4):
        if (r >> d) & 1: r ^= 0b10011 << (d - 4)
    return r
a5 = 1
for _ in range(5): a5 = gmul(a5, 2)
cols = [gmul(a5, 1 << j) for j in range(4)]         # images of basis vectors
A = tuple(tuple((cols[j] >> i) & 1 for j in range(4)) for i in range(4))
H = cyclic(A); print('alpha^5 =', a5, 'matrix rows', to_rows(A), '|H| =', len(H))
assert all(act(A, x) == gmul(a5, x) for x in DIRS)
T = [mul(h, tp(h)) for h in H]
for S in T: print('target', [list(r) for r in S])
want = [[[0,1,1,1],[1,0,1,1],[1,1,1,1],[1,1,1,0]], [[1,0,1,0],[0,1,1,0],[1,1,0,1],[0,0,1,1]]]
print('listed targets present:', all(tuple(map(tuple, w)) in T for w in want), 'I present:', ID in T)
def form(s):   # parse 'x1x2 + x1 + ...'
    Q = 0
    for t in s.replace(' ', '').split('+'):
        v = [int(c) - 1 for c in t.split('x') if c]
        i, j = (v[0], v[0]) if len(v) == 1 else (min(v), max(v))
        Q |= 1 << PAIRS.index((i, j))
    return Q
six = [('x1 + x1x4', [1,3,5,7]), ('x1x2 + x1x3 + x1x4 + x2x4 + x3x4 + x4', [3,5,8,14]), ('x2 + x2x3 + x2x4 + x3x4', [2,3,12,13]),
       ('x1x2 + x1x4 + x2x3 + x2x4 + x3x4 + x4', [3,6,8,13]), ('x1 + x1x2', [1,5,9,13]), ('x1x2 + x1x3 + x2 + x2x3', [2,5,10,13])]
g = Game(H, 4)
print('orbits', g.orbits)
for C, name in [(g.counts_from_dirs([1,2,4,8,3] * 3), 'pattern 1,2,4,8,3 x3'), (g.counts_from_dirs(DIRS), 'all 15 once')]:
    print(name, 'virtual moves per direction:', sorted(set(C.values())))
C = g.counts_from_dirs(DIRS)
sl = lambda Q: sum(C[x] for x in SUPP[Q]) - g.nodd[Q]
for s, sp in six:
    Q = form(s); print('%-45s supp %s (listed %s) odd at %d of 3 targets, slack %d' % (s, sorted(SUPP[Q]), sp, sum(pairing(Q, S) for S in T), sl(Q)))
ntight = [Q for Q in range(1, NF) if len(SUPP[Q]) == 4 and sum(pairing(Q, S) for S in T) == 3]
print('4-point forms odd at all 3 targets:', len(ntight), ' tight (slack 0) forms overall:', sum(1 for Q in range(1, NF) if sl(Q) == 0))
F = form('x1x2 + x1x3'); print('F supp', sorted(SUPP[F]), 'odd at', [pairing(F, S) for S in T], '|D_F| =', sum(C[x] for x in SUPP[F]), 'slack', sl(F))
Z = set()
for Q in range(1, NF):
    if sl(Q) == 0:
        for x in SUPP[Q]:
            for y in SUPP[Q]: Z.add((x, y))
print('blocks among 3,5,13 forced zero:', all((x, y) in Z for x in (3,5,13) for y in (3,5,13)))
print('blocks with 11:', [(x, y, (x, y) in Z) for x in (3,5,11,13) for y in (11,) ] + [(11, y, (11, y) in Z) for y in (3,5,13)])
print('verdict all-once:', g.verdict(C))
