#!/usr/bin/env python3
"""gridcover.py (w8-rmexact): rank cover with GRID partitions (independent row and column partitions of supp F) and
forced zeros, for fixed direction counts; supports of size <= --smax (default 6).  For a grid {R_i} x {C_j}:
slack(F) <= sum_{i,j} u(R_i, C_j), u = 0 if every cell of R_i x C_j is a forced zero, else
min(C(R_i), C(C_j), min_{G contains R_i u C_j} slack(G)).  Exact integers; numpy only for the sums.
Usage: gridcover.py --H singer:5 --comps "3,3,4,2,3" """
import sys, argparse, itertools, numpy as np
sys.path.insert(0, '.')
from forms import setup, pair
ap_ = argparse.ArgumentParser()
ap_.add_argument('--m', type=int, default=4); ap_.add_argument('--H', default='singer:5'); ap_.add_argument('--poly', type=int, default=None)
ap_.add_argument('--conj', type=int, default=-1); ap_.add_argument('--W', type=int, default=4)
ap_.add_argument('--comps', default='3,3,4,2,3'); ap_.add_argument('--smax', type=int, default=6)
args = ap_.parse_args(); m = args.m
if args.poly is None: args.poly = {3: 0b1011, 4: 0b10011, 5: 0b100101}[m]
D = setup(m); supp = D['supp']; nf = D['nf']; dirs = D['dirs']; mm, tr = D['matmul'], D['transpose']; ap = D['apply']
exec(open('hgroups.py').read())
H = make_H(args, D); n = len(H)
T = [D['enc_sym'](mm(h, tr(h))) for h in H]
nodd = [sum(pair(F, S) for S in T) for F in range(nf)]
orb = {}; orbs = []
for x in dirs:
    if x in orb: continue
    o = sorted({ap(h, x) for h in H})
    for y in o: orb[y] = len(orbs)
    orbs.append(o)
S = [frozenset(supp[F]) for F in range(nf)]
def set_partitions(lst):
    if not lst: yield []; return
    first, rest = lst[0], lst[1:]
    for p in set_partitions(rest):
        for i in range(len(p)): yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p
for comp in [tuple(int(v) for v in c.split(',')) for c in args.comps.split(';')]:
    C = {x: comp[orb[x]] for x in dirs}
    sl = [sum(C[x] for x in S[F]) - args.W * nodd[F] for F in range(nf)]
    if min(sl[1:]) < 0: print(comp, 'deletion'); continue
    tight = [F for F in range(1, nf) if sl[F] == 0 and nodd[F] > 0]
    zero = set()
    for F in tight:
        for x in S[F]:
            for y in S[F]: zero.add((x, y))
    msl = {}
    def ms(U):
        if U not in msl: msl[U] = min(sl[G] for G in range(1, nf) if U <= S[G])
        return msl[U]
    verdict = ('pass',)
    worst = None
    for F in range(1, nf):
        P = sorted(x for x in S[F] if C[x] > 0)
        if len(P) > args.smax or sl[F] == 0: continue
        subs = [frozenset(c) for r in range(1, len(P) + 1) for c in itertools.combinations(P, r)]
        si = {s_: i for i, s_ in enumerate(subs)}
        U = np.zeros((len(subs), len(subs)), dtype=np.int64)
        for a, R in enumerate(subs):
            for b, Cc in enumerate(subs):
                if all((x, y) in zero for x in R for y in Cc): U[a, b] = 0
                else: U[a, b] = min(sum(C[x] for x in R), sum(C[y] for y in Cc), ms(R | Cc))
        parts = list(set_partitions(P))
        Aidx = np.zeros((len(parts), len(subs)), dtype=np.int64)
        for k, p in enumerate(parts):
            for blk in p: Aidx[k, si[frozenset(blk)]] = 1
        best = None
        for k in range(len(parts)):
            V = Aidx[k] @ U                    # sum over row blocks of this row partition
            costs = Aidx @ V                   # every column partition
            c = int(costs.min())
            if best is None or c < best: best = c
        if best < sl[F]: verdict = ('grid-cover', F, P, 'needs', sl[F], 'cover', best); break
        if worst is None or best - sl[F] < worst[0]: worst = (best - sl[F], F, len(P))
    print(comp, verdict, 'tightest margin', worst, flush=True)
