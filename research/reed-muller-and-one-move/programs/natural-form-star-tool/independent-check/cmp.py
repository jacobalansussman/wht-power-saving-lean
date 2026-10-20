"""cmp.py <expect.txt> <out files...>: compare rm2-exact verdicts with vrec verdicts (letter and, for Z, round)."""
import sys
from collections import Counter
exp = {}
for l in open(sys.argv[1]):
    t = l.split(); exp[t[0]] = t[1]
c = Counter(); bad = []
for f in sys.argv[2:]:
    for l in open(f):
        t = l.split(' ', 1)
        if len(t) < 2 or t[0] not in exp: continue
        mine = t[1].strip()[2]; e = exp[t[0]]
        key = (e, mine); c[key] += 1
        ok = (e == mine) or (e in ('Z1', 'Z2') and mine == 'Z') or (e == '8a' and mine == 'C') or (e in ('P', 'deepX') and mine == 'P')
        if not ok: bad.append(l.strip())
print(dict(c), "disagreements:", len(bad)); [print(" ", b[:120]) for b in bad[:10]]
