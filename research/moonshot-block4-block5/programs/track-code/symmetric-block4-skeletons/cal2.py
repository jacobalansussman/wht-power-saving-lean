"""cal2.py -- calibration 2: the rotated plain schedule (array l uses directions e_{l}, e_{l+1}, e_{l+2}) is a union of
tau-orbits under Z3rot at m = 3 (W = 3) and under Z4rot at m = 4 (W = 4); check orbit sizes and the exact span test."""
import sys
from symsearch import *
for (m, W, Hn) in [(3, 3, 'Z3rot'), (4, 4, 'Z4rot'), (4, 2, 'Z2swap')]:
    R = Run(m, W, Hn, 1); G = R.G
    orders = [[1 << ((l + j) % m) for j in range(m)] for l in range(W)]
    gens = plain_program(R.S, W, orders)
    # orbit check: tau of every generator is (up to free) some generator
    ok = True
    for u in gens:
        for h in R.Hc:
            v = R.tau(h, u)
            if not any(R.S.same_up_to_free(v, w) for w in gens):
                ok = False
    print(f"m={m} W={W} H={Hn} |H|={len(R.H)}: plain rotated program tau-stable: {ok}; span test: {final_all(R.S, gens, W)}")
