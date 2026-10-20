"""openjoint2.py -- w5-twoshared: one more producer round on the Lemma-D-open cores.  For pairs (top-0 set, top-1 set)
of Lemma-D injections, the joint skeleton's closure cells are recomputed and every exact single producer and
producer pair (sweep2.producers) for them is offered; every subset of size <= 2 is decided (both outputs, every
coefficient free).  Usage: openjoint2.py seed core npairs maxpool cap"""
import sys, os, itertools, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
seed, cidx, npairs, maxpool, cap = [int(a) for a in sys.argv[1:6]]
sys.argv = ['opencase.py', str(seed), str(cidx), '1', str(cap)]
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'opencase.py')).read().split("# (i) each output alone")[0]
exec(compile(src, 'opencase_build', 'exec'))
import sweep2
from exact import Closure
tally = {}
pairs = list(itertools.product(range(len(inj_by_top[0])), range(len(inj_by_top[1]))))
random.Random(5).shuffle(pairs)
for (a, b) in pairs[:npairs]:
    base = list(inj_by_top[0][a][1]) + list(inj_by_top[1][b][1])
    mv0 = moves_of(R0, R1, None, None, base)
    sweep2.EXIST.clear()
    for mi, summ in enumerate(mv0): sweep2.EXIST[mi] = [(s_, M) for (s_, M) in summ]
    P0 = build(2, mv0)
    cells = []; sz = []
    for l0 in range(2):
        C = Closure(P0, l0, cap)
        if not C.ok: sz.append('cap'); continue
        sz.append((C.nr, C.nc)); inv_ = {v: k for k, v in C.rows.items()}
        cells += [inv_[r] for r in range(C.nr)]
    pool = sweep2.producers(P0, R0, R1, cells)[:maxpool]
    def expand(items):
        ex_ = []
        for it in items:
            if it[0] == 'pair': ex_ += [it[1], it[2]]
            else: ex_.append(it)
        return ex_
    res = {}
    for sub in [()] + [(i,) for i in range(len(pool))] + list(itertools.combinations(range(len(pool)), 2)):
        mvs = moves_of(R0, R1, None, None, base + expand([pool[i] for i in sub]))
        P = build(2, mvs)
        scr, sizes = singular_script(P, cap=cap)
        v = 'CAP' if scr is None else decide(scr, f'open_sing/j2_s{seed}_c{cidx}.sing', timeout=200)[0]
        res[v] = res.get(v, 0) + 1; tally[v] = tally.get(v, 0) + 1
        if v == 'CONSISTENT':
            pickle.dump({'moves': mvs}, open(f'open_sing/J2_CONS_s{seed}_c{cidx}_{a}_{b}_{"_".join(map(str, sub))}.pkl', 'wb'))
            print("   CONSISTENT subset", sub, sizes, flush=True)
    print(f"pair ({a},{b}): joint closures {sz}, producers {len(pool)}: {res}", flush=True)
print("TALLY", tally, flush=True)
