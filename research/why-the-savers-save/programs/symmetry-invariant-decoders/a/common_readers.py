"""common_readers.py (t-km, part a): for two face planes P (block b) and Q (block b' != b), the targets that know BOTH
(read P and read Q).  A register whose frame contains P + Q can serve only such targets (frame inside u_t^perp), so
this caps any cross-block sharing of one register between the factories of two blocks (projection argument).
Also the rank of P + Q.  usage: python3 -I common_readers.py"""
import json, itertools
from collections import Counter, defaultdict
d = json.load(open('a0_demand.json')); ports = d['ports']; v = len(ports)
dem = {int(t): set(int(s) for s in dd) for t, dd in d['demand'].items()}
blocks = json.load(open('u8_blocks.json'))['blocks']; blk = {s: b for b, f in enumerate(blocks) for s in f}
def rank(xs):
    b = []
    for x in xs:
        for y in b: x = min(x, x ^ y)
        if x: b.append(x); b.sort(reverse=True)
    return len(b)
readers = defaultdict(set)
for t in range(v):
    per = defaultdict(list)
    for s in dem[t]: per[blk[s]].append(s)
    for b, ss in per.items():
        if t not in blocks[b]: readers[(b, tuple(sorted(ss)))].add(t)
faces = sorted(readers)
hist = Counter(); rk = Counter()
for (f1, f2) in itertools.combinations(faces, 2):
    if f1[0] == f2[0]: continue
    c = len(readers[f1] & readers[f2])
    hist[c] += 1
    if c: rk[(c, rank([ports[s] for s in f1[1] + f2[1]]))] += 1
print('face pairs in different blocks by number of common readers:', sorted(hist.items()))
print('(common readers, rank of P+Q):', sorted(rk.items()))
per_face = Counter(sum(1 for g in faces if g[0] != f[0] and len(readers[f] & readers[g]) >= 4) for f in faces)
print('faces by number of partner faces with >= 4 common readers:', sorted(per_face.items()))
json.dump(dict(hist=sorted(hist.items()), rank=[[list(k), c] for k, c in sorted(rk.items())], partners4=sorted(per_face.items())),
          open('common_readers.json', 'w'))
