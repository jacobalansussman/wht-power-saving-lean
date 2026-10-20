#!/usr/bin/env python3
"""rcover.py (w8-rmexact): rectangle-cover rank bound with forced zeros, exact, for fixed direction counts.
For a form F with support P (|P| <= --smax), rank E[D_F, D_F] = slack(F) (Theorem 10').  E vanishes on cells (x,y)
inside a common tight support.  Cover the remaining cells of P x P by rectangles R x C (R, C subsets of P; a cell
is in exactly one rectangle; forced-zero cells may be in several or none); rank is subadditive over this sum, and rank E[R, C] <= u(R,C) = min(C(R), C(C),
min slack(G) over forms G with R u C in supp G).  So slack(F) <= min-cost cover.  Exact DP over subsets of cells.
Usage: rcover.py --H singer:5 --comps file|all [--smax 4]."""
import sys, argparse, itertools
sys.path.insert(0, '.')
from forms import setup, pair
ap_ = argparse.ArgumentParser()
ap_.add_argument('--m', type=int, default=4); ap_.add_argument('--H', default='singer:5'); ap_.add_argument('--poly', type=int, default=None)
ap_.add_argument('--conj', type=int, default=-1); ap_.add_argument('--W', type=int, default=4); ap_.add_argument('--M', type=int, default=15)
ap_.add_argument('--comps', default='all'); ap_.add_argument('--smax', type=int, default=4)
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
k = len(orbs)
S = [frozenset(supp[F]) for F in range(nf)]
def test(comp):
    C = {x: comp[orb[x]] for x in dirs}
    sl = [sum(C[x] for x in S[F]) - args.W * nodd[F] for F in range(nf)]
    if min(sl[1:]) < 0: return ('deletion',)
    tight = [F for F in range(1, nf) if sl[F] == 0 and nodd[F] > 0]
    zero = set()
    for F in tight:
        for x in S[F]:
            for y in S[F]: zero.add((x, y))
    msl = {}
    def ms(U):
        if U not in msl: msl[U] = min(sl[G] for G in range(1, nf) if U <= S[G])
        return msl[U]
    for F in range(1, nf):
        P = sorted(x for x in S[F] if C[x] > 0)
        if len(P) > args.smax or sl[F] == 0: continue
        cells = [(x, y) for x in P for y in P if (x, y) not in zero]
        if not cells:
            if sl[F] > 0: return ('cover', F, P, 'all cells zero', sl[F])
            continue
        ci = {c: i for i, c in enumerate(cells)}; full = (1 << len(cells)) - 1
        rects = []
        subsets = [frozenset(c) for r in range(1, len(P) + 1) for c in itertools.combinations(P, r)]
        for R in subsets:
            for Cc in subsets:
                mask = 0
                for x in R:
                    for y in Cc:
                        if (x, y) in ci: mask |= 1 << ci[(x, y)]
                if not mask: continue
                u = min(sum(C[x] for x in R), sum(C[y] for y in Cc), ms(R | Cc))
                rects.append((mask, u))
        # keep, per mask, the cheapest rectangle
        best = {}
        for mk, u in rects:
            if mk not in best or u < best[mk]: best[mk] = u
        rects = list(best.items())
        INF = 10 ** 9; dp = [INF] * (full + 1); dp[0] = 0
        for st in range(full + 1):
            if dp[st] == INF: continue
            # cover the lowest uncovered cell
            if st == full: break
            low = (~st) & (st + 1)
            for mk, u in rects:
                if mk & low and not (mk & st):          # a PARTITION of the nonzero cells (overlaps only on forced zeros)
                    ns = st | mk
                    if dp[st] + u < dp[ns]: dp[ns] = dp[st] + u
        if dp[full] < sl[F]: return ('rect-cover', F, P, 'needs', sl[F], 'cover', dp[full])
    return ('pass', len(tight))
if args.comps == 'all':
    comps = []
    for cut in itertools.combinations(range(args.M + k - 1), k - 1):
        b = (-1,) + cut + (args.M + k - 1,)
        comps.append(tuple(b[i + 1] - b[i] - 1 for i in range(k)))
else: comps = [tuple(int(v) for v in s.split(',')) for s in args.comps.split(';')]
npass = 0
for c in comps:
    r = test(c); npass += r[0] == 'pass'
    print(c, r, flush=True)
print("pass: %d of %d" % (npass, len(comps)))
