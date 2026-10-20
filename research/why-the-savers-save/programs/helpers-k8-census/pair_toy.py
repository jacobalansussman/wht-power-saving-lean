#!/usr/bin/env python3
"""g-break pair toy (H1 + H2 break): one Hamming block in abstract coordinates.
Block span S = F2^4; labels A = the 8 vectors with bit 3 set (an affine coset: checked on U8, probes.json P1); the
14 hyperplanes other than H0 = {bit3 = 0} hold 4 labels each; U8 uses 12 (one parallel pair unused).
Delivery item for a plane P: either a PACKET (row on P's 4 labels, created at frame P: k visitors at P) or TWO PAIR ROWS
(P's labels split into 2 pairs {a,b},{c,d}; each pair row lives on the line <a,b>, k visitors per line; a line serves
every plane that uses it).  k = 2 (one copy) or 3 (two copies), as in the hub's lemma, applied per frame.
Registers have nested chains: at most one line and one plane each, line inside plane if both (the hub's rule, per rank).
Sources: label in every frame they visit.  Helpers: anything.
Exact for a given design: min registers = units - max matching(line units, plane units containing them); helpers upper
bound = that cover minus a max matching of sources into its chains; lower bound = min registers - 8.
Search: exhaustive over the line set when packets are fixed, random restarts otherwise.  Usage: python3 -I pair_toy.py"""
import itertools, random, json, sys
from collections import defaultdict
A = [x for x in range(16) if x >> 3 & 1]
dot = lambda a, b: bin(a & b).count('1') & 1
hyper = {f: frozenset(x for x in A if not dot(x, f)) for f in range(1, 16) if f != 8}
drop = {1, 9}            # one parallel pair (f, f+8) unused, as in U8 (12 planes)
planes = [f for f in hyper if f not in drop]
lines = [frozenset(p) for p in itertools.combinations(A, 2)]
splits = {f: [] for f in planes}
for f in planes:
    S = sorted(hyper[f])
    for a in S[1:]:
        L1 = frozenset((S[0], a)); L2 = frozenset(set(S) - L1); splits[f].append((L1, L2))
def matching(left, adj):
    match = {}
    def aug(u, seen):
        for w in adj[u]:
            if w in seen: continue
            seen.add(w)
            if w not in match or aug(match[w], seen): match[w] = u; return True
        return False
    return sum(1 for u in left if aug(u, set())), match
def evaluate(choice, k):
    """choice: plane -> None (packet) or split index. returns (lower, upper) helpers"""
    L = set(); Pk = []
    for f, ch in choice.items():
        if ch is None: Pk.append(f)
        else: L |= set(splits[f][ch])
    lu = [(l, i) for l in L for i in range(k)]; pu = [(f, i) for f in Pk for i in range(k)]
    adj = {u: [p for p in pu if u[0] <= hyper[p[0]]] for u in lu}
    m, match = matching(lu, adj)
    regs = len(lu) + len(pu) - m
    chains = [(match[p][0], p[0]) for p in match] + [(u[0], None) for u in lu if u not in match.values()] + \
             [(None, p[0]) for p in pu if p not in match]
    sadj = {s: [ci for ci, (l, f) in enumerate(chains) if (l is None or s in l) and (f is None or s in hyper[f])] for s in A}
    ms, _ = matching(A, sadj)
    return max(regs - 8, 0), regs - ms, len(L), len(Pk)
res = {}
random.seed(1)
for k in (2, 3):
    allpk = evaluate({f: None for f in planes}, k)
    best = None
    for trial in range(4000):
        ch = {f: (None if random.random() < random.random() else random.randrange(3)) for f in planes}
        # local improvement
        cur = evaluate(ch, k)
        improved = True
        while improved:
            improved = False
            for f in planes:
                for opt in [None, 0, 1, 2]:
                    if opt == ch[f]: continue
                    old = ch[f]; ch[f] = opt; e = evaluate(ch, k)
                    if (e[1], e[0]) < (cur[1], cur[0]): cur = e; improved = True
                    else: ch[f] = old
        if best is None or (cur[1], cur[0]) < (best[0][1], best[0][0]):
            best = (cur, dict(ch))
    (lo, up, nl, npk), ch = best
    res[f'k={k}'] = dict(all_packets_helpers=allpk[1], best_upper_helpers=up, its_lower=lo, lines=nl, packet_planes=npk,
                         design={str(f): ch[f] for f in planes})
    print(k, res[f'k={k}'])
json.dump(res, open('pair_toy.json', 'w'), indent=1)
