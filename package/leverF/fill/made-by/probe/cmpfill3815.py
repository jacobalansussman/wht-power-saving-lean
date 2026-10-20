# o-model: compare two staged544 T360 owner fills (our code). Widths from PR3734's NATIVE-base (as #3815 build_fill.py).
# Reports: literal heads (empty bins of children 8/10), remainders per child-pair class, number of classes whose
# multiset of bin patterns differs, and number of bins whose (pattern) differs.  Args: NATIVE-base.json.gz TEMPLATE FILL_A FILL_B
import json, gzip, sys
from collections import Counter, defaultdict
rd = lambda p: json.loads(gzip.decompress(open(p, 'rb').read())) if p.endswith('.gz') else json.load(open(p))
def rank(U):
    piv = {}
    for x in U:
        while x:
            b = x.bit_length()
            if b in piv: x ^= piv[b]
            else: piv[b] = x; break
    return len(piv)
nat = rd(sys.argv[1]); tmpl = rd(sys.argv[2])
width = [rank(nat['frames'][nat['final'][r]]) - rank(nat['frames'][nat['start'][r]]) for r in range(224, 632)]
caps = {c['child']: c['rank'] for c in tmpl['idle_children']}
pairs = [(0, 2), (1, 3), (4, 6), (5, 7), (8, 10), (9, 11), (12, 13), (14, 15), (16, 17)]
pair_of = {c: i for i, cs in enumerate(pairs) for c in cs}
def load(p):
    rows = rd(p)['rows']; bins = defaultdict(list)
    for t, h, o, c, pt, z in rows: bins[o, c, pt].append(width[h])
    return {k: tuple(sorted(v)) for k, v in bins.items()}, rows
A, ra = load(sys.argv[3]); B, rb = load(sys.argv[4])
T = 360
def stats(X):
    heads = sum(1 for o in range(T) for c in (8, 10) for p in range(112) if (o, c, p) not in X)
    rem = defaultdict(Counter)
    for o in range(T):
        for c, cap in caps.items():
            for p in range(112):
                rem['%d/%d' % pairs[pair_of[c]]][cap - sum(X.get((o, c, p), ()))] += 1
    return heads, {k: dict(sorted(v.items())) for k, v in sorted(rem.items())}
for name, X in (('A', A), ('B', B)):
    h, rem = stats(X); print(name, 'heads', h); print('  remainders by pair class', rem)
cls_diff = 0
for o in range(T):
    for cs in pairs:
        ma = Counter(A.get((o, c, p), ()) for c in cs for p in range(112)); mb = Counter(B.get((o, c, p), ()) for c in cs for p in range(112))
        cls_diff += ma != mb
print('classes (owner, child pair) with a different multiset of bin patterns:', cls_diff, 'of', 9 * T)
print('bins with a different pattern:', sum(1 for o in range(T) for c in caps for p in range(112) if A.get((o, c, p), ()) != B.get((o, c, p), ())))
print('rows with same (t,h,o,c,p):', sum(1 for x, y in zip(ra, rb) if x[:5] == y[:5]), 'of', len(ra))
print('rows with same (t,h) -> same (o, child pair):', sum(1 for x, y in zip(ra, rb) if x[2] == y[2] and pair_of[x[3]] == pair_of[y[3]]))
