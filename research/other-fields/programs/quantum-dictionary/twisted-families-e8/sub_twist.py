#!/usr/bin/env python3
"""sub_twist.py (t-family) -- uniform-twist cost of the coordinate (frame) decoder on sub-systems of the E8 lines:
E7 (lines orthogonal to 2e_0), E8 minus axes (the hub's P112), D8 (the 56 lines e_i +- e_j, scaled), and D4.
Per total a: hub labels' cost = F_2-rank of B0 = (1 + G) mod 2 on A_a (exists only where B0 restricted is
the right pattern: B0 = 1 exactly at orthogonal pairs, true for E8 lines); twisted lower bound = clique of pairwise
non-orthogonal lines in A_a; exact twisted minimum = least k with a completion (CSP of e8_minrank7.py through a
k-clique), tried for k = clique only.  python3 -I sub_twist.py"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import e8_minrank7 as E
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '../lib')))
from copyfloor import rank
L, G = E.L, E.G
out = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)
def run(name, lines, coords):
    tot_old = tot_lo = tot_tw = 0; v = len(lines)
    for a in coords:
        A = [S for S in lines if L[S][a] != 0]
        if len(A) < 2: continue
        adj = {S: {T for T in A if T != S and G[S][T] != 0} for S in A}
        cl = E.max_clique(A, adj); c = len(cl)
        B0 = [[(1 + G[i][j]) % 2 for j in A] for i in A]
        old = rank(B0, 2)
        Gs = E.G
        r, n, dt, ds, asg = E.solve(A, cl, adj)
        tw = c if r else '>%d' % c
        tot_old += old; tot_lo += c; tot_tw += c if r else c + 1
        say('  %s total %d: |A_a| %d, clique %d, hub-label rank %d, twisted completion of rank %d: %s'
            % (name, a, len(A), c, old, c, 'yes' if r else 'no'))
    say('%s: v %d; hub-label cst %d (x %.4f); twisted cst >= %d, best found %s (x %.4f)'
        % (name, v, tot_old, tot_old / v, tot_lo, tot_tw, tot_tw / v))
e7 = [i for i in range(120) if G[i][0] == 0]
run('E7 (63 lines orthogonal to 2e_0)', e7, range(1, 8))
run('E8 minus axes (P112)', list(range(8, 120)), range(8))
# D4 inside E8: lines with support in coordinates {0,1,2,3} -> the 1,1,1,1-type lines and axes
d4 = [i for i in range(120) if all(L[i][c] == 0 for c in range(4, 8))]
run('E8 lines inside coordinates 0..3 (%d lines)' % len(d4), d4, range(4))
open(os.path.join(HERE, 'out', 'sub_twist.txt'), 'w').write('\n'.join(out) + '\n')
