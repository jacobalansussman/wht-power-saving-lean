"""sweep_cap.py -- w5-twoshared (resumed): decide the skeletons that sweep2.py / sweep_sig.py left as CAP (closure over
60 columns).  Replays the same cores (same seeds, same rng sequence) and the same producer pools, then
 (1) decides the UNION of the whole pool (every coefficient free: INCONSISTENT for the union decides every subset,
     because a subset's solution extends by zeros), if its closures fit in maxcol columns;
 (2) otherwise decides each subset that was CAP at cap 60 OUTPUT BY OUTPUT at cap2 (one output INCONSISTENT excludes
     the skeleton; a unit needs both outputs).
Prints the per-core CAP count at cap 60 so the replay can be compared with the original output.
Usage: sweep_cap.py variant ncores seed maxpool cap2 maxcol timeout outdir     (variant: sig1, sig0, cross, priv)"""
import sys, os, random, itertools, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sweep2 import *
import sweep2
from sweep_sig import core_sig

variant, ncores, seed, maxpool, cap2, maxcol, tmo, outdir = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), sys.argv[8]
os.makedirs(outdir, exist_ok=True)
rng = random.Random(seed); tally = {}
CAP0 = 60
for c in range(ncores):
    t0 = time.time()
    if variant in ('sig1', 'sig0'):
        cr = core_sig(rng, variant)
        if cr is None: print("core failed", flush=True); continue
        R0, R1, a11, a20, tg, o = cr
    else:
        cr = core(rng, variant)
        if cr is None: print("core failed", flush=True); continue
        R0, R1, a11, a20, tg = cr
    mv0 = moves_of(R0, R1, a11, a20)
    sweep2.EXIST.clear()
    for mi, summ in enumerate(mv0):
        sweep2.EXIST[mi] = [(src, Mx) for (src, Mx) in summ]
    P0 = build(2, mv0)
    cells = []; sz = []
    for l0 in range(2):
        C = Closure(P0, l0, CAP0)
        if not C.ok: sz.append('cap'); continue
        sz.append((C.nr, C.nc))
        inv = {v: k for k, v in C.rows.items()}
        cells += [inv[r] for r in range(C.nr)]
    pool = producers(P0, R0, R1, cells)
    npool = len(pool)
    pool = pool[:maxpool]
    def expand(items):
        ex = []
        for it in items:
            if it[0] == 'pair': ex += [it[1], it[2]]
            else: ex.append(it)
        return ex
    subsets = [()] + [(i,) for i in range(len(pool))] + list(itertools.combinations(range(len(pool)), 2))
    capsubs = []
    for sub in subsets:
        P = build(2, moves_of(R0, R1, a11, a20, expand([pool[i] for i in sub])))
        if not all(Closure(P, l, CAP0).ok for l in range(2)):
            capsubs.append(sub)
    print(f"core {c}: tags {tg}, bare closures {sz}, producers {npool}, pool {len(pool)}, {len(subsets)} skeletons, CAP at {CAP0}: {len(capsubs)}", flush=True)
    if not capsubs:
        continue
    # (1) the union of the whole pool, output by output (a solution of a subset, restricted to one output, extends by
    #     zeros to the union; so output l of the union INCONSISTENT => output l of every subset INCONSISTENT)
    def per_output(mvs, tag):
        P = build(2, mvs); info = []
        for l in range(2):
            C = Closure(P, l, cap2)
            if not C.ok or C.nc > maxcol:
                info.append((l, 'cap' if not C.ok else (C.nr, C.nc), 'CAP')); continue
            scr, sizes = singular_script(P, cap=cap2, outputs=[l])
            v, out = decide(scr, os.path.join(outdir, f'{tag}.sing'), timeout=tmo)
            info.append((l, sizes[0], v))
            if v == 'INCONSISTENT': break
        return info
    info = per_output(moves_of(R0, R1, a11, a20, expand(pool)), f'U_{variant}_{seed}_{c}')
    print(f"   union of the pool, per output: {info} ({time.time()-t0:.0f}s)", flush=True)
    if any(v == 'INCONSISTENT' for (_, _, v) in info):
        tally['INCONSISTENT'] = tally.get('INCONSISTENT', 0) + len(capsubs)
        print(f"core {c}: all {len(capsubs)} CAP skeletons INCONSISTENT (via one output of the union)", flush=True)
        continue
    # (2) each CAP subset, output by output
    res = {}
    for sub in capsubs:
        info = per_output(moves_of(R0, R1, a11, a20, expand([pool[i] for i in sub])), f'{variant}_{seed}_{c}')
        v = 'INCONSISTENT' if any(x[2] == 'INCONSISTENT' for x in info) else ('CONSISTENT-PARTIAL' if any(x[2] == 'CONSISTENT' for x in info) else 'CAP')
        res[v] = res.get(v, 0) + 1; tally[v] = tally.get(v, 0) + 1
        print(f"   subset {sub}: {info} -> {v} ({time.time()-t0:.0f}s)", flush=True)
    print(f"core {c}: CAP skeletons re-decided: {res} ({time.time()-t0:.0f}s)", flush=True)
print("TOTAL", tally, flush=True)
