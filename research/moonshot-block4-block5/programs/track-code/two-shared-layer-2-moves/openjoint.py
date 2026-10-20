"""openjoint.py -- w5-twoshared: JOINT decision for the Lemma-D-open private cores.  "Output e alone" is not a test
(with every coefficient free the solver may drop the other array from v); here both outputs are decided together on
skeletons carrying BOTH tops' Lemma-D injection sets (top 0's three and top 1's three), every coefficient free.
Usage: openjoint.py seed core nlift n0 n1 cap"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib, pickle
seed, cidx, nlift, n0, n1, cap = [int(a) for a in sys.argv[1:7]]
sys.argv = ['opencase.py', str(seed), str(cidx), str(nlift), str(cap), 'nodecide']
# reuse opencase's construction without its decisions
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'opencase.py')).read()
src = src.split("# (i) each output alone")[0]
exec(compile(src, 'opencase_build', 'exec'))
tally = {}
import itertools
for a, b in itertools.product(range(min(n0, len(inj_by_top[0]))), range(min(n1, len(inj_by_top[1])))):
    tt0, new0 = inj_by_top[0][a]; tt1, new1 = inj_by_top[1][b]
    mvs = moves_of(R0, R1, None, None, list(new0) + list(new1))
    P = build(2, mvs)
    scr, sizes = singular_script(P, cap=cap)
    if scr is None: v = 'CAP'
    else: v, out = decide(scr, f'open_sing/joint_s{seed}_c{cidx}.sing', timeout=300)
    tally[v] = tally.get(v, 0) + 1
    print(f"joint: top0 set {a} {tt0}, top1 set {b} {tt1}: closures {sizes} {v}", flush=True)
    if v == 'CONSISTENT':
        pickle.dump({'moves': mvs}, open(f'open_sing/JOINT_CONS_s{seed}_c{cidx}_{a}_{b}.pkl', 'wb'))
        print(out[-2000:], flush=True)
print("TALLY", tally, flush=True)
