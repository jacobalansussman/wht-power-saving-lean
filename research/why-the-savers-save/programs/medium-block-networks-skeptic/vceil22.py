#!/usr/bin/env python3
"""v-cx-medium / vceil22.py (own code; python3 -I).  Generous ceilings for a one-pair design at block m = 1 + 3 sigma,
in cx-medium's own most generous model (junk cleaning free, data pieces {sigma, 2 sigma}, centre copies as t pieces of
rank c/t per call), with my fold (vprice.fold_root), and the same with U8X0's helper mass PER JUNK ENTRY on ideal banks
(U8X0: 10,030 narrow per call at m = 25 for 5,824 junk entries; scaled by ln(m/r)/ln(25/r) per piece)."""
import sys, math, json
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from vprice import fold_root, cls
def ceiling(n, sigma, c, t, junk_per_port):
    m = 1 + 3 * sigma
    H = {sigma: 2 * n, 2 * sigma: 2 * n}
    r = c / t
    H[r] = H.get(r, 0) + 3 * t
    D = m * 2 * n - sum(k * v for k, v in H.items())
    a_free = fold_root(H, 2 * n, m) if D > 0 else 0.0
    # helpers at U8X0's per-entry piece mix
    per_entry = {rr: nn / 5824 for rr, nn in cls['s'].items()}
    Hh = dict(H)
    jun = n * junk_per_port
    lev = 0
    for rr, f in per_entry.items():
        Hh[rr] = Hh.get(rr, 0) + 3 * jun * f; lev += 3 * jun * f * rr
    W = 2 * n + lev / m
    a_h = fold_root(Hh, W, m) if D > 0 else 0.0
    return dict(m=m, D=D, a_junk_free=a_free, a_U8_helper_density=a_h)
if __name__ == '__main__':
    n, sigma, c, t, jp = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    print(json.dumps(ceiling(n, sigma, c, t, jp)))
