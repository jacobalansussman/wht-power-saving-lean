"""witness.py (d-relay) -- for a saved no-closure relay solution, print relayed pairs (S, T) whose relay frame G
cannot exist: T's register is clean (y_T + x_T) only after T cancelled x_T' at a frame containing a_T', and a_T'
is not orthogonal to a_S, so no frame inside a_S^perp follows that cancellation on T's chain.
usage: python3 -I witness.py out/d5_L3_none.json"""
import sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import patterns
d = json.load(open(sys.argv[1]))
knows, G = patterns.instance(d['instance'])
Ls = patterns.lines_of([[1 if k in (i, j) else 0 for k in range(5)] for i in range(5) for j in range(5) if i < j] +
                       [[1 if k == i else -1 if k == j else 0 for k in range(5)] for i in range(5) for j in range(5) if i < j]) if d['instance'] == 'D5' else None
if d['instance'] == 'D5':
    import itertools
    rs = []
    for i, j in itertools.combinations(range(5), 2):
        for s in (1, -1):
            v = [0] * 5; v[i] = 1; v[j] = s; rs.append(v)
    Ls = patterns.lines_of(rs)
X, Y, lay = d['X'], d['Y'], d['layer']; v = len(X); r = len(X[0])
B = [[sum(X[S][k] & Y[T][k] for k in range(r)) & 1 for T in range(v)] for S in range(v)]
print('instance', d['instance'], 'rank', r, 'layers', lay)
n = 0
for S in range(v):
    for T in range(v):
        if S == T or knows[S][T] or not B[S][T]: continue
        bad = [Tp for Tp in range(v) if Tp != T and knows[T][Tp] and B[T][Tp] and not knows[S][Tp]]
        if bad and n < 3:
            Tp = bad[0]; n += 1
            print('S=%d %s relays T=%d %s (layer %d < %d); T cancels T\'=%d %s (a_T\'.a_T = %d) but a_T\'.a_S = %d != 0'
                  % (S, Ls[S] if Ls else '', T, Ls[T] if Ls else '', lay[T], lay[S], Tp, Ls[Tp] if Ls else '', G[T][Tp], G[S][Tp]))
print('relayed pairs whose frame G cannot exist (one-level test):',
      sum(1 for S in range(v) for T in range(v) if S != T and not knows[S][T] and B[S][T] and any(Tp != T and knows[T][Tp] and B[T][Tp] and not knows[S][Tp] for Tp in range(v))),
      'of', sum(1 for S in range(v) for T in range(v) if S != T and not knows[S][T] and B[S][T]))
