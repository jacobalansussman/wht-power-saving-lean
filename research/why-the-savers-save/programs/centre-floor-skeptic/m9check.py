#!/usr/bin/env python3
"""v-centre / m9check.py: own check of the m = 9 (profile 7^9) and m = 10 (1,6,7^8) cover claims.
Co-flats of size <= 72: Gamma_a (112, 53), N_i (8, 56), D7-type with no F0 line (128, 70), rebuilt here."""
import json, itertools, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'data.json')))
lines = [tuple(v) for v in D['lines']]
Gam = [int(x) for x in D['Gam']]; Ni = [int(x) for x in D['Ni']]
code = set()
G = [0b11110000, 0b11001100, 0b10101010, 0b11111111]
for c in itertools.product([0, 1], repeat=4):
    w = 0
    for b, g in zip(c, G):
        if b: w ^= g
    code.add(w)
F0 = [tuple(2 if j == i else 0 for j in range(8)) for i in range(8)]
ip = lambda a, b: sum(x * y for x, y in zip(a, b))
d7 = set()
for x in itertools.product(range(-2, 3), repeat=8):
    if sum(v * v for v in x) != 8: continue
    if sum((v % 2) << (7 - i) for i, v in enumerate(x)) not in code: continue
    cof = sum(1 << s for s in range(112) if ip(lines[s], x) != 0)
    if bin(cof).count('1') <= 72: d7.add(cof)
print('D7 co-flats <= 72 ports:', len(d7), set(bin(c).count('1') for c in d7))
cofl = Gam + Ni + sorted(d7)
L7 = sorted(set(int(l, 16) for l in open(os.path.join(HERE, 'L7.txt'))))
L6 = None
pc = lambda x: bin(x).count('1')
print('dim-7 closed sets:', len(L7), 'sizes', sorted(set(pc(s) for s in L7)))
edges = set()
for C in cofl:
    c = pc(C)
    # partner covers at most max_s |s & C|
    inter = [(pc(s & C), i) for i, s in enumerate(L7)]
    mx = max(k for k, _ in inter)
    cand = [i for k, i in inter if k >= c - mx]
    for a, b in itertools.combinations(cand, 2):
        if (L7[a] | L7[b]) & C == C: edges.add((a, b))
print('covering-pair edges:', len(edges))
adj = {}
for a, b in edges: adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
tri = sum(1 for a, b in edges for c in adj[a] & adj[b] if c > b)
print('triangles in the covering-pair graph:', tri, '(a 9-clique needs triangles)')
# m = 10 profile (1,6,7^8): a port + a dim-6 set + a dim-7 set covering Gamma_a or N_i (57 < 70 rules out others)
L6 = [int(l, 16) for l in open(os.path.join(HERE, 'L6.txt'))]
L6 = sorted(set(L6))
hits = 0
for C in Gam + Ni:
    c = pc(C)
    c7 = [s for s in L7 if pc(s & C) >= c - 1 - 20]
    c6 = [s for s in L6 if pc(s & C) >= c - 1 - 32]
    for s7 in c7:
        rem = C & ~s7
        for s6 in c6:
            if pc(rem & ~s6) <= 1: hits += 1
print('m=10 (1,6,7^8) port+dim6+dim7 covers of Gamma_a/N_i:', hits)
