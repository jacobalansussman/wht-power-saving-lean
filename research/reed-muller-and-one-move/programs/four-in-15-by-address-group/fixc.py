#!/usr/bin/env python3
"""fixc.py (w8-rmexact): the multi-target rank-cover test with the virtual direction counts FIXED by the physical
directions.  With c fixed the LP of mtlp.py decouples: it is feasible iff (i) slack(G) = c(G) - a_G >= 0 for every G,
and (ii) for every form F and every partition {R_i} of supp F used,
      sum_{i,j} u(R_i, R_j) >= slack(F),   u(R,C) = min(c(R), c(C), min_{G: R u C in supp G} slack(G)).
(Setting every piece variable to its upper bound is optimal.)  Exact arithmetic (fractions).
H-orbit compositions: for H acting freely on directions (Singer powers), c_x = n_o / (W |H|) with n_o the number of
physical moves in the orbit o of x.  Usage: fixc.py --H singer:5 --W 4 --M 15 [--allpart 4] [--all | --comp 3,3,3,3,3]
Prints every composition (up to nothing) with its verdict and the first violated row."""
import sys, argparse, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from forms import setup, pair
ap_ = argparse.ArgumentParser()
ap_.add_argument('--m', type=int, default=4); ap_.add_argument('--H', default='singer:5')
ap_.add_argument('--poly', type=int, default=None); ap_.add_argument('--conj', type=int, default=-1)
ap_.add_argument('--W', type=int, default=4); ap_.add_argument('--M', type=int, default=15)
ap_.add_argument('--allpart', type=int, default=4); ap_.add_argument('--all', action='store_true')
ap_.add_argument('--comp', default='')
ap_.add_argument('--part', default='0/1')
args = ap_.parse_args(); m = args.m
if args.poly is None: args.poly = {3: 0b1011, 4: 0b10011, 5: 0b100101}[m]
D = setup(m); supp = D['supp']; nf = D['nf']; dirs = D['dirs']; mm, tr = D['matmul'], D['transpose']
exec(open('hgroups.py').read())
H = make_H(args, D); n = len(H)
T = [D['enc_sym'](mm(h, tr(h))) for h in H]
aF = [Fr(sum(pair(F, S) for S in T), n) for F in range(nf)]
# direction orbits under H
ap = D['apply']; orb = {}; orbs = []
for x in dirs:
    if x in orb: continue
    o = sorted({ap(h, x) for h in H}); 
    for y in o: orb[y] = len(orbs)
    orbs.append(o)
assert all(len(o) == n for o in orbs), "H must act freely on directions"
S = [frozenset(supp[F]) for F in range(nf)]
def set_partitions(lst):
    if not lst: yield []; return
    first, rest = lst[0], lst[1:]
    for p in set_partitions(rest):
        for i in range(len(p)): yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p
parts = {}
for F in range(1, nf):
    P = sorted(S[F])
    parts[F] = list(set_partitions(P)) if len(P) <= args.allpart else [[[x] for x in P]]
# all unions R u C that occur
unions = set()
for F, pl in parts.items():
    for part in pl:
        for X in part:
            for Y in part: unions.add(frozenset(X) | frozenset(Y))
contain = {U: [G for G in range(1, nf) if U <= S[G]] for U in unions}
def check(comp):
    c = {x: Fr(comp[orb[x]], args.W * n) for x in dirs}
    sl = [sum(c[x] for x in S[G]) - aF[G] for G in range(nf)]
    for G in range(1, nf):
        if sl[G] < 0: return ('deletion', G, sl[G])
    ms = {U: min(sl[G] for G in contain[U]) for U in unions}
    worst = None
    for F in range(1, nf):
        for part in parts[F]:
            tot = Fr(0)
            for X in part:
                cX = sum(c[x] for x in X)
                for Y in part:
                    cY = sum(c[y] for y in Y)
                    tot += min(cX, cY, ms[frozenset(X) | frozenset(Y)])
            gap = tot - sl[F]
            if gap < 0: return ('cover', F, [len(X) for X in part], gap)
            if worst is None or gap < worst[0]: worst = (gap, F, [len(X) for X in part])
    return ('pass', 'min margin', worst)
k = len(orbs)
comps = []
if args.comp: comps = [tuple(int(v) for v in args.comp.split(','))]
elif args.all:
    for cut in itertools.combinations(range(args.M + k - 1), k - 1):
        b = (-1,) + cut + (args.M + k - 1,)
        comps.append(tuple(b[i + 1] - b[i] - 1 for i in range(k)))
pi_, pn_ = map(int, args.part.split('/')); comps = comps[pi_::pn_]
print("H=%s conj=%d |H|=%d orbits of directions %s; W=%d M=%d allpart=%d; %d compositions" % (args.H, args.conj, n, [o[0] for o in orbs], args.W, args.M, args.allpart, len(comps)), flush=True)
npass = 0
for comp in comps:
    r = check(comp)
    if r[0] == 'pass': npass += 1
    if r[0] == 'pass' or len(comps) <= 20: print(comp, r, flush=True)
print("compositions passing: %d of %d" % (npass, len(comps)))
