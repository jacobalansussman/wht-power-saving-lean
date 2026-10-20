# independent exact replay of a word with arbitrary declared start frames (no start-at-0 rule); frames as label sets,
# nested + strict dim increase on every climb, F2 symbolic values; final checks as xcheck_p Part F.
import gzip, json, sys
from fractions import Fraction as Q
from collections import Counter
w = json.load(gzip.open(sys.argv[1], 'rt')); w0 = json.load(gzip.open(sys.argv[2], 'rt'))
V, R = w['v'], w['R']; N = 2*V + R; L = w['lines']
def rk(vs):
    B = []
    for v in vs:
        r = [Q(x) for x in v]
        for p, b in B:
            if r[p]: c = r[p]; r = [x - c*y for x, y in zip(r, b)]
        p = next((i for i, x in enumerate(r) if x), None)
        if p is not None: c = r[p]; B.append((p, [x / c for x in r]))
    return len(B)
fs = [frozenset(f) for f in w['frames']]; fd = [rk([L[i] for i in f]) for f in w['frames']]
cur = list(w['start']); val = [1 << i for i in range(N)]; H = Counter()
def climb(r, f):
    if cur[r] == f: return
    assert fs[cur[r]] <= fs[f] and fd[cur[r]] < fd[f], ('bad climb', r, cur[r], f)
    H[fd[f] - fd[cur[r]]] += 1; cur[r] = f
def run(sec):
    for f, t, ss in w[sec]:
        for r in [t, *ss]: climb(r, f)
        for s in ss: val[t] ^= val[s]
run('D'); run('A')
tops = [fd[cur[q]] for q in w['ret']]; cp = [val[q] for q in w['ret']]
for t, row in enumerate(w['scatter']):
    for k in row: val[V+t] ^= cp[k]
run('B'); run('C')
for r in range(N): climb(r, w['final'][r])
ok = all(val[s] == 1 << s for s in range(V)) and all(val[h] == 1 << h for h in range(2*V, N)) and all(val[V+t] == (1 << (V+t)) | (1 << t) for t in range(V))
paid = H + Counter(tops)
print('values exact:', ok, 'paid', dict(sorted(paid.items())), 'mass', sum(k*n for k, n in paid.items()))
print('D reads 444:', any(444 in ss for f, t, ss in w['D']), '| start dims of helpers', Counter(fd[w['start'][h]] for h in range(2*V, N)))
print('final helper dims', Counter(fd[w['final'][h]] for h in range(2*V, N)))
dd = [(i, g0, g) for i, (g0, g) in enumerate(zip(w0['D'], w['D'])) if g0 != g]
print('D gates changed:', len(dd))
ia = [i for i in range(min(len(w0['A']), len(w['A']))) if w0['A'][i] != w['A'][i]][:1]
if ia: i = ia[0]; print('first A diff at', i, w0['A'][i], '->', w['A'][i])
