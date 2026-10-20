"""Find reads (gate sources) whose row is held, at that moment, by another register too (up to a dyadic scalar)."""
import sys, gzip, json, importlib.util
sys.path.insert(0,'.')
from fractions import Fraction as F
from collections import Counter, defaultdict
spec = importlib.util.spec_from_file_location('al', 'audit_local_py39.py'); al = importlib.util.module_from_spec(spec); spec.loader.exec_module(al)
import fw
w = fw.load('K8-FULL-WORD.json.gz'); v=112; R=480; n=2*v+R
ident = lambda: [{k: F(1)} for k in range(n)]
main = w['A'] + w['cut'] + w['B']
raw = ident()
for e in main: al.apply(e, raw)
dirt = [{} for _ in range(n)]
for y in range(v, 2*v): dirt[y] = {k: a for k, a in raw[y].items() if k >= 2*v}
for e in w['precondition']: al.apply(e, dirt)
comp = [['add', y, k, -a.numerator, a.denominator, []] for y in range(v, 2*v) for k, a in sorted(dirt[y].items())]
# rows as linear forms in the initial values: row[r] = dict col->coef. Applying gate t += q s : row[t] += q row[s]
rows = ident()
def key(d):
    if not d: return ('zero',)
    # normalise by first coefficient (projective) to catch dyadic multiples
    k0 = min(d); c = d[k0]
    return tuple(sorted((k, x / c) for k, x in d.items())), 
def full_key(d): return tuple(sorted(d.items()))
for e in w['precondition'] + comp: al.apply(e, rows)
hk = {r: key(rows[r]) for r in range(n)}
groups = defaultdict(set)
for r in range(n): groups[hk[r]].add(r)
stat = Counter(); ex = []
secs = [('A', e) for e in w['A']] + [('cut', e) for e in w['cut']] + [('B', e) for e in w['B']]
for idx, (sec, e) in enumerate(secs):
    if e[0] == 'add' and sec != 'cut':
        s = e[2]; others = groups[hk[s]] - {s, e[1]}
        if others:
            stat[(sec, 'read-has-twin')] += 1
            if len(ex) < 12: ex.append((sec, idx, e[:5], sorted(others)[:6]))
        else: stat[(sec, 'read-unique')] += 1
    t = e[1]
    old = hk[t]
    al.apply(e, rows)
    groups[old].discard(t); hk[t] = key(rows[t]); groups[hk[t]].add(t)
print(stat)
for x in ex: print(x)
