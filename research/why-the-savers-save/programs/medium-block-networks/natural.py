#!/usr/bin/env python3
"""cx-medium / natural.py (own code, stdlib; python3 -I).  Costs of the coordinate (natural) decoders: an orthonormal
basis of the design's span, totals T_f = sum_b (f . r_b) x_b; per total: support size, hub cost dim span(u_S) and
pairing rank rank B[S,S].  E8/P112: e_1..e_8; E7 (span = (1,1,0..)^perp): (1,-1,0,..), e_3..e_8; D_n: e_i;
A_n: the n+1 ambient e_i (a tight frame of the sum-zero hyperplane: sum_i (e_i.r)(e_i.s) = <r,s>)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from designs_cx import design, labels, ip, rank2
from fractions import Fraction
out = {}
for name in ('P112', 'E8', 'E7', 'E6', 'D8', 'D6', 'A5'):
    L = design(name); labs, k, B = labels(L); n = len(L); dim = len(L[0]); norm = ip(L[0], L[0])
    if name == 'E7':
        F = [tuple([1, -1] + [0] * 6)] + [tuple(1 if j == i else 0 for j in range(8)) for i in range(2, 8)]
    elif name == 'E6':
        F = [tuple([1, 1, 1] + [0] * 5)] + [tuple(1 if j == i else 0 for j in range(8)) for i in range(3, 8)]
    else:
        F = [tuple(1 if j == i else 0 for j in range(dim)) for i in range(dim)]
    # check the delivered map: M = c * sum_f (f.r_a)(f.r_b) must be delta on B = 1 pairs
    G = [[ip(f, l) for l in L] for f in F]
    M = [[sum(Fraction(G[q][a] * G[q][b], 1) / Fraction(ip(F[q], F[q])) for q in range(len(F))) for b in range(n)] for a in range(n)]
    c = 1 / M[0][0]
    ok = all((c * M[a][b] == (1 if a == b else 0)) for a in range(n) for b in range(n) if B[a][b] == 1)
    junk = [abs(c * M[a][b]) for a in range(n) for b in range(n) if B[a][b] == 0 and M[a][b] != 0]
    rows = []
    for q, f in enumerate(F):
        S = [a for a in range(n) if G[q][a] != 0]
        hub = rank2([labs[a] for a in S])
        pr = rank2([int(''.join(str(B[a][b]) for b in S), 2) for a in S])
        rows.append(dict(f=f, support=len(S), hub=hub, pair=pr))
    hubc = sum(r['hub'] for r in rows); prc = sum(r['pair'] for r in rows)
    out[name] = dict(n=n, k=k, totals=len(F), exact_on_B1=ok, junk_entries=len(junk), junk_values=sorted(set(map(str, junk))),
                     rows=rows, hub_cost=hubc, pair_cost=prc, be_P1=2 * n / 3, be_P2=4 * n / 5)
    print(f"{name}: n={n} k={k} totals={len(F)} exact on B=1 pairs {ok}; junk entries {len(junk)} values {sorted(set(map(str, junk)))};"
          f" cost per call hub {hubc} ({hubc/n:.3f}/port) pair {prc} ({prc/n:.3f}/port); break-even P=1 {2*n/3:.1f} P=2 {4*n/5:.1f};"
          f" per total (support, hub, pair): {[(r['support'], r['hub'], r['pair']) for r in rows]}")
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'natural.json'), 'w'), indent=1, default=str)
