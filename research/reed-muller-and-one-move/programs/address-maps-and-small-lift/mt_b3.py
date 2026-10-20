"""mt_b3.py (rm-theory): block 3 probe of target diversity. Multi-target slice = no-phase Singer game (Theorem B).
Shapes: 5 in 14 (2.8 = the hook bound, not excluded by any test in the canonical placement), controls 5 in 15 plain
(must be found), and the impossible 2 in 5 (Theorem W2) and 3 in 8 (below the hook 2.8) which must NOT be found."""
import sys, time
import numpy as np
from mt_lm import Setup, MTSystem
from mt_cover import cover_test

def runinst(label, m, Wspec, dspec, order, restarts, iters, seed):
    st = Setup(m, Wspec, dspec, order)
    et = st.exact_tests(); ct = cover_test(st)
    if et["excluded"] or ct["excluded"]:
        print(f"{label}: {Wspec} {dspec} {order}: EXCLUDED exactly {et} {ct}", flush=True); return None
    sy = MTSystem(st); rng = np.random.default_rng(seed); best = 1e99; hits = 0
    for r in range(restarts):
        x0 = (rng.standard_normal(sy.nvar) + 1j * rng.standard_normal(sy.nvar)) / np.sqrt(2); x0[2 * sy.nA:] *= 0.3
        x, f, it = sy.solve(x0, iters=iters)
        best = min(best, f); hits += f < 1e-18
    print(f"{label}: {Wspec} {dspec} {order}: W={st.W} M={st.M} tight={et['tight']} best residual {best:.3e}, hits {hits}/{restarts}", flush=True)
    return best

def negatives(restarts, iters):
    for ts, comp in [("0,4,5", (1, 1, 2, 3, 4, 5, 6, 7)), ("0,1,3", (1, 2, 2, 3, 4, 5, 6, 7)), ("0,2,6", (1, 2, 3, 4, 4, 5, 6, 7)),
                     ("0,0,1,3", (1, 1, 1, 2, 2, 3, 4, 4, 5, 6, 7)), ("0,4,5,5", (1, 1, 1, 2, 3, 3, 4, 5, 6, 6, 7)),
                     ("0,4,4,5", (1, 1, 1, 2, 3, 4, 5, 5, 6, 7, 7))]:
        for order in ["lex", "rand1"]:
            runinst("IMPOSSIBLE (hook 2.8)", 3, "k:" + ts, "list:" + ",".join(map(str, comp)), order, restarts, iters, 7)


if __name__ == "__main__":
    restarts = int(sys.argv[1]); iters = int(sys.argv[2]); which = sys.argv[3]
    d14 = 'list:' + ','.join(str(z) for z in list(range(1, 8)) * 2)
    if which == "neg":
        negatives(restarts, iters)
    elif which == "controls":
        # plain 5 in 15: frames of 5 targets
        runinst("control plain 5 in 15", 3, "k:0,1,2,3,4", "frames", "lex", restarts, iters, 1)
        runinst("control plain 5 in 15 rand", 3, "k:0,1,2,3,4", "frames", "rand3", restarts, iters, 2)
        for ks in ["k:0,1", "k:0,3", "k:0,0"]:
            for dl in ["list:1,2,3,4,5", "list:1,2,4,7,3", "list:1,2,3,5,6", "list:1,1,2,4,7"]:
                runinst("impossible 2 in 5", 3, ks, dl, "lex", restarts, iters, 3)
        for ks in ["k:0,1,2", "k:0,1,3"]:
            for dl in ["list:1,2,3,4,5,6,7,1", "list:1,2,4,3,5,6,7,7"]:
                runinst("impossible 3 in 8", 3, ks, dl, "lex", restarts, iters, 4)
    else:
        for ks in ["k:0,1,2,3,4", "k:0,1,2,4,5", "k:0,1,3,4,5", "k:0,0,1,2,4", "k:0,1,2,3,5"]:
            for order in ["lex", "singer", "rand1", "rand2"]:
                runinst("5 in 14 (U_3^(2))", 3, ks, d14, order, restarts, iters, 5)

