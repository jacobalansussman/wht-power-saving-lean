# root cut-LP dual for each pi (from pis file); writes "pi... | y y_1 ... y_254" with y scaled to be exactly dual feasible
import sys
sys.path.insert(0, '.')
import numpy as np, cutip
from scipy.optimize import linprog
src, dst = sys.argv[1], sys.argv[2]
out = open(dst, 'w')
for line in open(src):
    v = line.split('#')[0].split()
    if len(v) < 8: continue
    pi = tuple(int(x) for x in v[:8])
    r = cutip.rvec(pi)
    res = linprog(np.ones(len(cutip.arcs)), A_ub=-cutip.A, b_ub=-r, bounds=(0, None), method="highs")
    y = np.maximum(-res.ineqlin.marginals, 0)
    col = (y @ cutip.A).max()
    if col > 1: y = y / col
    y = np.floor(y * 1e9) / 1e9          # round down: stays feasible
    out.write(' '.join(map(str, pi)) + ' | y ' + ' '.join(f'{t:.9f}' for t in y) + f'   # LP {res.fun:.6f} dualbound {float(y @ r):.6f}\n')
out.close()
