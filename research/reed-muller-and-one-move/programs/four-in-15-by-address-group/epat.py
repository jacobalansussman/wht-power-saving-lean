#!/usr/bin/env python3
"""epat.py (w8-rmexact): ex2-families Corollary 10' (E-pattern matching test) on the virtual multi-target system.
E[i,j] = 0 is forced when virtual moves i, j (i = j allowed) lie in a common tight support (Theorem 10'(a): on a tight
support E[D,D] = 0).  Theorem 10': rank E[D_F, D_F] = |D_F| - n_odd(F) for every F, so the bipartite graph
(rows D_F) x (cols D_F) minus the forced zeros must have a matching of that size.  The zero pattern depends only on
directions, so the matching is a b-matching on directions with capacities C_x (max flow).  Exact, integer.
Usage: epat.py --H singer:5 --W 4 --comps <file with one composition per line, or 'all'>"""
import sys, argparse, itertools
sys.path.insert(0, '.')
from forms import setup, pair
ap_ = argparse.ArgumentParser()
ap_.add_argument('--m', type=int, default=4); ap_.add_argument('--H', default='singer:5'); ap_.add_argument('--poly', type=int, default=None)
ap_.add_argument('--conj', type=int, default=-1); ap_.add_argument('--W', type=int, default=4); ap_.add_argument('--M', type=int, default=15)
ap_.add_argument('--comps', default='all')
ap_.add_argument('--dirs', default='', help='explicit physical directions (overrides --comps)')
ap_.add_argument('--verbose', action='store_true')
args = ap_.parse_args(); m = args.m
if args.poly is None: args.poly = {3: 0b1011, 4: 0b10011, 5: 0b100101}[m]
D = setup(m); supp = D['supp']; nf = D['nf']; dirs = D['dirs']; mm, tr = D['matmul'], D['transpose']; ap = D['apply']
exec(open('hgroups.py').read())
H = make_H(args, D); n = len(H)
T = [D['enc_sym'](mm(h, tr(h))) for h in H]
nodd = [sum(pair(F, S) for S in T) for F in range(nf)]          # odd channels; odd virtual arrays = W * nodd
orb = {}; orbs = []
for x in dirs:
    if x in orb: continue
    o = sorted({ap(h, x) for h in H})
    for y in o: orb[y] = len(orbs)
    orbs.append(o)
k = len(orbs)
def maxflow(cap, s, t, N):
    flow = 0
    while True:
        par = [-1] * N; par[s] = s; st = [s]
        while st and par[t] < 0:
            u = st.pop()
            for v in range(N):
                if par[v] < 0 and cap[u][v] > 0: par[v] = u; st.append(v)
        if par[t] < 0: return flow
        v = t; b = float('inf')
        while v != s: u = par[v]; b = min(b, cap[u][v]); v = u
        v = t
        while v != s: u = par[v]; cap[u][v] -= b; cap[v][u] += b; v = u
        flow += b
def test(comp):
    C = {x: comp[orb[x]] for x in dirs}      # virtual moves per direction (free action: n_o per direction)
    sl = [sum(C[x] for x in supp[F]) - args.W * nodd[F] for F in range(nf)]
    if min(sl[1:]) < 0: return ('deletion',)
    tightF = [F for F in range(1, nf) if sl[F] == 0 and nodd[F] > 0]
    zero = set()
    for F in tightF:
        P = sorted(supp[F])
        for x in P:
            for y in P: zero.add((x, y))
    worst = None
    for F in range(1, nf):
        P = sorted(x for x in supp[F] if C[x] > 0)
        need = sl[F]
        if need == 0: continue
        N = 2 * len(P) + 2; s, t = N - 2, N - 1
        cap = [[0] * N for _ in range(N)]
        for i, x in enumerate(P):
            cap[s][i] = C[x]; cap[len(P) + i][t] = C[x]
            for j, y in enumerate(P):
                if (x, y) not in zero: cap[i][len(P) + j] = 10 ** 9
        f = maxflow(cap, s, t, N)
        if f < need: return ('E-pattern', F, sorted(supp[F]), 'needs rank', need, 'max matching', f)
        if worst is None or f - need < worst[0]: worst = (f - need, F)
    return ('pass', len(tightF), 'tight forms; min matching margin', worst)
comps = []
if args.dirs:
    dl = [int(v) for v in args.dirs.split(',')]; c = [0] * k
    for z in dl: c[orb[z]] += 1
    comps = [tuple(c)]; args.M = len(dl)
elif args.comps == 'all':
    for cut in itertools.combinations(range(args.M + k - 1), k - 1):
        b = (-1,) + cut + (args.M + k - 1,)
        comps.append(tuple(b[i + 1] - b[i] - 1 for i in range(k)))
else:
    comps = [tuple(int(v) for v in line.split(',')) for line in open(args.comps) if line.strip()]
print("H=%s |H|=%d W=%d M=%d: %d compositions" % (args.H, n, args.W, args.M, len(comps)), flush=True)
npass = 0
for c in comps:
    r = test(c)
    if r[0] == 'pass': npass += 1; print(c, r, flush=True)
    elif args.verbose or len(comps) <= 5: print(c, r, flush=True)
print("pass deletion + E-pattern: %d of %d" % (npass, len(comps)))
