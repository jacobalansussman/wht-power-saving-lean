"""exact_sweep.py -- w3-asym: exact (Groebner) sweep of small skeletons around Lemma-5 cores.
For each core ((2,2) single-letter core, or the sigma-merge core of famS): a pool of injections = the landing letters
(x_j into u read through v_0 / v_1, both outputs) + frame-matched injections for the cells of the bare core's closures
(famI.candidates).  Every subset of size 1 and 2 of the pool (and size 3 of the landing letters) is built, its closures
computed (cap), and its polynomial system (program + readout coefficients, both outputs) decided by Singular std().
Usage: exact_sweep.py kind ncores seed nsamp cap outdir   (kind: 22 or sigma)"""
import sys, os, itertools, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_decide import *
from famS import core_sigma

def decide(scr, fn):
    open(fn, 'w').write(scr)
    try:
        out = subprocess.run(['Singular', '-q', fn], capture_output=True, text=True, timeout=40).stdout
    except subprocess.TimeoutExpired:
        return 'TIMEOUT'
    if 'INCONSISTENT' in out: return 'INCONSISTENT'
    if 'CONSISTENT' in out: return 'CONSISTENT ' + out.strip().split('\n')[-1]
    return 'ERR'

if __name__ == '__main__':
    kind, ncores, seed, nsamp, cap, outdir = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)
    W = 2; tally = {}
    for s in range(ncores):
        rec = core22(rng, W) if kind == '22' else core_sigma(rng, W)
        if rec is None: print("core failed", flush=True); continue
        P0 = build_rec(W, rec)
        u = 2 * W
        pm = {}
        for (l, g) in P0.terms[u]:
            if G.dist(g) == 2 and l not in pm: pm[l] = G.MATS[g]
        ells = [rec[W + 1][0][1], rec[W + 3][0][1]]
        pool = []
        for (l, e) in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            j = 1 - l
            J = G.mul(G.K, np.conj(pm[l]).T, pm[j]); fk = G.frame_key(J)
            for k in landing_letters(rng, fk, G.mul(G.D1, ells[e]), np.eye(16, dtype=complex), nsamp):
                pool.append((W, j, k, f'land x{j}->u out{l} v{e}'))
        nland = len(pool)
        Cs = [Closure(P0, l0, cap) for l0 in range(W)]
        cells = []
        for C in Cs:
            if C.ok:
                cells += bad_cells(C, np.ones(P0.np, dtype=complex), top=3)
        cands = []
        for cell in cells[:2]:
            cands += candidates(P0, rec, W, cell, rng, 60)
        rng.shuffle(cands)
        for (mi, src, b, how) in cands[:6]:
            pool.append((mi, src, b, f'cand move{mi} src g{src} {how}'))
        subsets = [()] + [(i,) for i in range(len(pool))] + list(itertools.combinations(range(len(pool)), 2)) + list(itertools.combinations(range(nland), 3))
        res_core = {}
        for sub in subsets:
            r2 = [list(m) for m in rec]
            for i in sub:
                mi, src, b, name = pool[i]
                r2[mi].append((src, b))
            P = build_rec(W, r2)
            scr, sizes = singular_script(P, cap=cap)
            if scr is None:
                verdict = 'CAP'
            else:
                verdict = decide(scr, os.path.join(outdir, f'{kind}_{seed}_{s}.sing'))
            res_core[verdict.split()[0]] = res_core.get(verdict.split()[0], 0) + 1
            tally[verdict.split()[0]] = tally.get(verdict.split()[0], 0) + 1
            if not verdict.startswith('INCONS') and verdict != 'CAP':
                print(f"   core {s} subset {[pool[i][3] for i in sub]} sizes {sizes}: {verdict}", flush=True)
        print(f"core {s} ({kind}): pool {len(pool)} ({nland} landing), {len(subsets)} skeletons: {res_core}", flush=True)
    print("TOTAL", tally, flush=True)
