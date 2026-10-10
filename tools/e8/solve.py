"""solve.py -- the solver stage: choose a whole hierarchy from a pool of sums with the integer programme of xmilp.py,
with two device weights in its objective.

The programme minimises  uses + ordinary sums.  An ordinary sum is formed on the helper of one part while the other
part's helper is left unused.  Two kinds of ordinary sum are counted less than 1, because a later pass can remove
the helper they cost:
   lam   an ordinary sum with a single-source part (the source array itself can play the part: xsub.py);
   mu    an ordinary sum whose smaller part has 2 .. MUK sources (xexp.py can replace such a helper by direct reads).
The certificate of this repository was made with lam = 0.55, mu = 0.85, no symmetry imposed (group 'trivial'),
time limit 30 s, tie-break seed 0 (= costs not disturbed).

Needs numpy and scipy: scipy.optimize.milp calls the HiGHS solver that is shipped inside scipy.  The recorded run
used Python 3.9.6, numpy 2.0.2, scipy 1.13.1.  Another solver version may return another optimum of equal value,
and then the passes after it give another certificate (see README.md).

usage: solve.py <pool.json> <lam> <mu> <out hierarchy.json> [time limit s = 30] [tie-break seed = 0] [group = trivial]"""
import sys, time, random
sys.dont_write_bytecode = True
import xmilp
from xlib import *
MUK = 2

def main(poolf, lam, mu, outf, tlim=30.0, seed=0, gname="trivial"):
    T0 = time.time()
    A = Act(gname)
    pool = hio.load_pool(poolf)
    M = xmilp.build(pool, A, log=lambda *a: None)
    syms, D, Q, so_size = M["syms"], M["D"], M["Q"], M["so_size"]
    def nsrc(i): return bin(syms[i][0]).count("1")
    print("model built: sums %d, group %s  %.0fs" % (len(syms) - V, gname, time.time() - T0)); sys.stdout.flush()
    t = time.time()
    def extra(cost, oX, oA, oY, oC, mu=mu, seed=seed):
        if seed:
            rr = random.Random(seed)
            for i in range(len(cost)):
                if cost[i] > 0: cost[i] += 1e-3 * rr.random()
        if mu >= 1.0: return
        for q in range(len(Q)):
            ds = [d for d in Q[q][1] if D[d][4] == 1]
            for d in Q[q][1]:
                w, u, v = D[d][5]
                small = min(nsrc(u), nsrc(v))
                if D[d][4] > 1 and 2 <= small <= MUK: cost[oA + d] = so_size[D[d][1]] * mu
            if not ds: continue
            w, u, v = D[ds[0]][5]
            small = min(nsrc(u), nsrc(v))
            if not (2 <= small <= MUK): continue
            wq = Q[q][0] * mu
            if len(ds) == 1:
                if D[ds[0]][0] != 2 * Q[q][0]: cost[oA + ds[0]] = wq
            else:
                cost[oC + q] = wq
    sol, info = xmilp.solve(M, lam, tlim, log=print, extra_cost=extra)
    if sol is None:
        print("lam %s mu %s: no solution" % (lam, mu)); return 1
    h = xmilp.hierarchy(M, sol, A, clean_first=True)
    pr = xmilp.proxy(h)
    xmilp.save_h(h, outf)
    print("seed %d lam %s mu %s: objective %.1f gap %s solve %.0fs | deliveries %d ordinary %d with source %d pairs %d | %s  %.0fs" % (
        seed, lam, mu, info["fun"], info.get("gap"), time.time() - t, pr["deliveries"], pr["ordinary"], pr["ordinary_with_source"], pr["pairs"], outf, time.time() - T0))
    return 0

if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) < 4:
        print(__doc__); sys.exit(2)
    sys.exit(main(a[0], float(a[1]), float(a[2]), a[3], float(a[4]) if len(a) > 4 else 30.0, int(a[5]) if len(a) > 5 else 0, a[6] if len(a) > 6 else "trivial"))
