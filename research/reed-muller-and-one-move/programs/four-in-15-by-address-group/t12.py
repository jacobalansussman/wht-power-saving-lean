#!/usr/bin/env python3
"""t12.py (w8-rmexact): ex2-families Theorem 12 (slack-s covering, via Lemma 11) on the virtual multi-target system,
fixed direction counts.  If a form Q_e with n_odd(Q_e) = 0 (E[D,D] nonsingular, Theorem 10') has |D| >= s+1 and
every set of at most s+2 of its virtual moves lies inside the support of a form G with slack(G) <= s (so every
principal submatrix of order <= s+2 has rank <= s), Lemma 11 makes E[D,D] singular: contradiction.  A set of virtual
moves lies in supp G iff its direction set does; a direction set U of size u needs u <= s+2 moves.
Also the odd version (Theorem 10'(c) with slack): not used.  Usage: t12.py --H singer:5 --comps "3,3,4,2,3" --smax 4"""
import sys, argparse, itertools
sys.path.insert(0, '.')
from forms import setup, pair
ap_ = argparse.ArgumentParser()
ap_.add_argument('--m', type=int, default=4); ap_.add_argument('--H', default='singer:5'); ap_.add_argument('--poly', type=int, default=None)
ap_.add_argument('--conj', type=int, default=-1); ap_.add_argument('--W', type=int, default=4)
ap_.add_argument('--comps', default='3,3,4,2,3'); ap_.add_argument('--smax', type=int, default=4)
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
for comp in [tuple(int(v) for v in c.split(',')) for c in args.comps.split(';')]:
    C = {x: comp[orb[x]] for x in dirs}
    sl = [sum(C[x] for x in S[F]) - args.W * nodd[F] for F in range(nf)]
    even = [F for F in range(1, nf) if nodd[F] == 0 and sum(C[x] for x in S[F]) > 0]
    print(comp, "min slack over odd forms:", min(sl[F] for F in range(1, nf) if nodd[F] > 0), "; forms even at every target:", len(even))
    found = None
    for s in range(0, args.smax + 1):
        good = [F for F in range(1, nf) if sl[F] <= s]
        for Fe in even:
            P = sorted(x for x in S[Fe] if C[x] > 0); Dn = sum(C[x] for x in P)
            if Dn < s + 1: continue
            ok = True
            for u in range(1, min(s + 2, len(P)) + 1):
                for U in itertools.combinations(P, u):
                    # need at least u moves, at most s+2: always possible since every direction has >= 1 move
                    Us = frozenset(U)
                    if not any(Us <= S[G] for G in good): ok = False; break
                if not ok: break
            if ok: found = (s, Fe, P); break
        if found: break
    print("   Theorem 12 verdict:", ("EXCLUDED with s = %d by the even form %d on directions %s" % found) if found else "no exclusion up to s = %d" % args.smax)
