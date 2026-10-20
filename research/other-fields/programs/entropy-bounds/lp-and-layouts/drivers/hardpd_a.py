# hardpd.py -- session 1's DIRECTED cut LP (PDLP with a rigorous dual bound, layout_cut_pdlp.py) on the two
# hill-climbed block-4 layouts whose three-end images fall short (10.75 and 11.00).
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # modules of ../
import time, dpg, layout_cut_pdlp as L
m, W = 4, 3
verts, idx, adj, dist = dpg.build(m); n = len(verts); arcs = [(u, v) for u in range(n) for v in adj[u]]
t0 = time.time()
for D in ([(331, 751), (1319, 596), (837, 1071)], [(2030, 1063), (1397, 1085), (1249, 145)]):
    A, c, lo, hi, nvar, npat = L.build(n, arcs, D, W * m)
    primal, LB, st = L.solve(A, c, lo, hi, nvar, W * m, time_limit=3000)
    print(f"{D}: PDLP {st}, primal {primal:.5f}, rigorous LB {LB:.5f} (routing 12) ({time.time()-t0:.0f}s)", flush=True)
