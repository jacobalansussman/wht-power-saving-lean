"""v-w5-twoshared pools.py: replay w5-twoshared's sweep cores (their builders) and print each core's producer pool
TYPES (move, source) truncated to maxpool = 8, to see whether skeletons use injections outside the single-letter model
of the frame tests: (4, 4|5) = u_k -> v (two-term gamma_k), (2|3, 2|3) = s_j -> u_k (two-term alpha_kj)."""
import sys, os, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'their'))
import sweep2
from sweep2 import *
from sweep_sig import core_sig
variant = sys.argv[1]; ncores = int(sys.argv[2]); seeds = [int(x) for x in sys.argv[3:]]
for seed in seeds:
    rng = random.Random(seed)
    for c in range(ncores):
        cr = core(rng, variant) if variant == 'cross' else core_sig(rng, variant)
        if cr is None: print('core failed', flush=True); continue
        R0, R1, a11, a20, tg = cr[:5]
        mv0 = moves_of(R0, R1, a11, a20)
        sweep2.EXIST.clear()
        for mi, summ in enumerate(mv0): sweep2.EXIST[mi] = [(src, Mx) for (src, Mx) in summ]
        P0 = build(2, mv0); cells = []
        for l0 in range(2):
            C = Closure(P0, l0, 60)
            if not C.ok: continue
            inv = {v: k for k, v in C.rows.items()}; cells += [inv[r] for r in range(C.nr)]
        pool = producers(P0, R0, R1, cells)[:8]
        typ = [('pair', it[1][0], it[1][1]) if it[0] == 'pair' else (it[0], it[1]) for it in pool]
        bad = [i for i, t in enumerate(typ) if t[0] != 'pair' and ((t[0] == 4 and t[1] in (4, 5)) or (t[0] in (2, 3) and t[1] in (2, 3)))]
        print(f"{variant} seed {seed} core {c}: pool types {typ}; outside single-letter model at indices {bad}", flush=True)
