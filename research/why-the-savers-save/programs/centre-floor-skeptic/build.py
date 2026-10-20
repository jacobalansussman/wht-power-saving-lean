#!/usr/bin/env python3
"""v-centre / build.py: own construction of P112 (Hamming coordinates), the dense pairing B, F2 labels, the orthogonal
8-frames, the E7-type co-flats and the D7-type hyperplanes.  Exports data.json for the C enumerator.  Own code."""
import itertools, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
# extended Hamming [8,4,4]: generator rows
G = [0b11110000, 0b11001100, 0b10101010, 0b11111111]
code = set()
for c in itertools.product([0, 1], repeat=4):
    w = 0
    for b, g in zip(c, G):
        if b: w ^= g
    code.add(w)
tet = sorted(w for w in code if bin(w).count('1') == 4)
assert len(tet) == 14
def bits(w): return [i for i in range(8) if w >> (7 - i) & 1]
lines = []
for w in tet:
    S = bits(w)
    for sg in itertools.product([1, -1], repeat=3):   # first nonzero coordinate +1: one representative per line
        v = [0] * 8; v[S[0]] = 1
        for i, s in zip(S[1:], sg): v[i] = s
        lines.append(tuple(v))
assert len(lines) == 112
F0 = [tuple(2 if j == i else 0 for j in range(8)) for i in range(8)]
def ip(a, b): return sum(x * y for x, y in zip(a, b))
# check E8: all 240 roots have norm 4, ips in {0,+-2,+-4}
allr = lines + F0
for a in allr:
    for b in allr:
        assert ip(a, b) in (0, 2, -2, 4, -4) and (ip(a, b) in (4, -4)) == (a == b)
n = 112
B = [[(1 if (t == s or ip(lines[t], lines[s]) == 0) else 0) for s in range(n)] for t in range(n)]
# check B_ts = (<r_t,r_s>/2 + 1) mod 2
for t in range(n):
    for s in range(n):
        assert B[t][s] == ((ip(lines[t], lines[s]) // 2 + 1) % 2)
rows = [sum(B[t][s] << s for s in range(n)) for t in range(n)]
# F2 row space basis by elimination, then coordinates of each row
basis = []   # list of (pivotbit, vec, combo) reduced
piv = {}
for r in rows:
    x = r
    for p, (vec) in sorted(piv.items(), reverse=True):
        if x >> p & 1: x ^= vec
    if x:
        p = x.bit_length() - 1
        # keep fully reduced
        for q in list(piv):
            if piv[q] >> p & 1: piv[q] ^= x
        piv[p] = x
rank = len(piv)
print('F2 rank of B:', rank)
assert rank == 9
pivs = sorted(piv)
def coords(r):
    x = r; c = 0
    for k, p in enumerate(pivs):
        if x >> p & 1:
            x ^= piv[p]; c |= 1 << k
    assert x == 0
    return c
lab = [coords(r) for r in rows]
assert len(set(lab)) == 112
# frames: 8 pairwise orthogonal lines
orth = [sum(1 << s for s in range(n) if s != t and ip(lines[t], lines[s]) == 0) for t in range(n)]
frames = []
def rec(cl, cand, start):
    if len(cl) == 8:
        frames.append(tuple(cl)); return
    c = cand
    while c:
        s = (c & -c).bit_length() - 1
        c &= c - 1
        if s < start: continue
        rec(cl + [s], cand & orth[s] & ~((1 << (s + 1)) - 1), s + 1)
for t in range(n):
    rec([t], orth[t] & ~((1 << (t + 1)) - 1), t + 1)
print('frames:', len(frames))
# co-flats of E7 type: Gamma_a = ports non-orthogonal to a (incl. a); N_i = ports with coordinate i nonzero
Gam = [sum(1 << s for s in range(n) if ip(lines[a], lines[s]) != 0) for a in range(n)]
Ni = [sum(1 << s for s in range(n) if lines[s][i] != 0) for i in range(8)]
print('Gamma sizes', set(bin(x).count('1') for x in Gam), 'N sizes', set(bin(x).count('1') for x in Ni))
# D7-type hyperplanes: normals = norm-8 lattice vectors x in Z^8, x mod 2 in code
d7 = []
for x in itertools.product(range(-2, 3), repeat=8):
    if sum(v * v for v in x) != 8: continue
    m2 = sum((v % 2) << (7 - i) for i, v in enumerate(x))
    if m2 not in code: continue
    d7.append(x)
print('norm-8 lattice vectors:', len(d7))
hist = {}
for x in d7:
    nr = sum(1 for r in allr if ip(r, x) == 0)
    f = sum(1 for r in F0 if ip(r, x) == 0)
    p = sum(1 for r in lines if ip(r, x) == 0)
    assert nr == 42 and p == 42 - f
    hist[112 - p] = hist.get(112 - p, 0) + 1
print('D7 co-flat sizes (vector count, each hyperplane twice):', dict(sorted(hist.items())))
json.dump(dict(lines=lines, lab=lab, frames=frames, Gam=[str(x) for x in Gam], Ni=[str(x) for x in Ni]),
          open(os.path.join(HERE, 'data.json'), 'w'))
with open(os.path.join(HERE, 'data.txt'), 'w') as f:
    f.write(' '.join(map(str, lab)) + '\n')
    f.write(str(len(frames)) + '\n')
    for fr in frames: f.write(' '.join(map(str, fr)) + '\n')
    cof = Gam + Ni
    f.write(str(len(cof)) + '\n')
    for c in cof: f.write(' '.join(str(s) for s in range(n) if c >> s & 1) + '\n')
