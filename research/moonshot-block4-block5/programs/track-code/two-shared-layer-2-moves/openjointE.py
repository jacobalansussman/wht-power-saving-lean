"""openjointE.py -- exact (Singular) joint decision for given pairs of Lemma-D injection sets whose joint closures are
small; also every pair (up to npairs) whose joint closures have <= maxcol columns per output.
Usage: openjointE.py seed core npairs maxcol cap"""
import sys, os, itertools, pickle, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
seed, cidx, npairs, maxcol, cap = [int(a) for a in sys.argv[1:6]]
sys.argv = ['opencase.py', str(seed), str(cidx), '1', str(cap)]
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'opencase.py')).read().split("# (i) each output alone")[0]
exec(compile(src, 'opencase_build', 'exec'))
from exact import Closure
pairs = list(itertools.product(range(len(inj_by_top[0])), range(len(inj_by_top[1]))))
random.Random(5).shuffle(pairs)
tally = {}; done = 0
for (a, b) in pairs:
    if done >= npairs: break
    mvs = moves_of(R0, R1, None, None, list(inj_by_top[0][a][1]) + list(inj_by_top[1][b][1]))
    P = build(2, mvs)
    Cs = [Closure(P, l, cap) for l in range(2)]
    if not all(C.ok for C in Cs) or max(C.nc for C in Cs) > maxcol: continue
    scr, sizes = singular_script(P, cap=cap)
    v, out = decide(scr, f'open_sing/jE_s{seed}_c{cidx}.sing', timeout=600)
    tally[v] = tally.get(v, 0) + 1; done += 1
    print(f"pair ({a},{b}) types {inj_by_top[0][a][0]} / {inj_by_top[1][b][0]}: closures {sizes} {v}", flush=True)
    if v == 'CONSISTENT':
        pickle.dump({'moves': mvs}, open(f'open_sing/JE_CONS_s{seed}_c{cidx}_{a}_{b}.pkl', 'wb')); print(out[-1500:], flush=True)
print("TALLY", tally, flush=True)
