"""ctl.py -- controls for exact.py in the two-shared builder: (a) plain 2 x 4 (8 moves) built as s, u (private), v_0, v_1
(private), t_0, t_1 -> CONSISTENT; (b) the same letters with v_0, v_1 merged into ONE shared v (7 moves: the private
variant of the 2+2+1+2 shape, plain routes) -> expected INCONSISTENT; (c) (a) with an extra free cross summand."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact import *
from asym import plain
P8 = plain(1)
L = [rec[0][1] for rec in P8.moves]          # the four plain letters (array 0); same for array 1
os.makedirs('sing_ctl', exist_ok=True)
# generators: 0 x0, 1 x1, 2 s0, 3 s1, 4 u1, 5 u2, ...
a = [(0, L[0])], [(1, L[0])], [(2, L[1])], [(3, L[1])]
mv8 = list(a) + [[(4, L[2])], [(5, L[2])], [(6, L[3])], [(7, L[3])]]
mv7 = list(a) + [[(4, L[2]), (5, L[2])], [(6, L[3])], [(6, L[3])]]
mv8x = list(a) + [[(4, L[2]), (5, L[2])], [(5, L[2])], [(6, L[3])], [(7, L[3])]]
for name, mv in (('plain8', mv8), ('merged7', mv7), ('plain8cross', mv8x)):
    P = build(2, mv)
    scr, sizes = singular_script(P)
    if scr is None: print(name, 'closure cap'); continue
    v, out = decide(scr, f'sing_ctl/{name}.sing')
    print(name, 'params', P.np, 'closures', sizes, v, flush=True)
