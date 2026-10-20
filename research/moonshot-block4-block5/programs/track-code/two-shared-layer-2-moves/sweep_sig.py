"""sweep_sig.py -- w5-twoshared: the SIGMA-PAIR variant of the 2+2+1+2 shape (BRIEF5: "a sigma-pair of shared moves,
then one shared layer-3 move").  u_1 = D(a10 s0 + a11 s1), u_2 = D(a10 s0 + a21 s1) (same array-0 letter), v = D(g1 u1 + g2 u2)
with g2 = D^-1 sig D g1, sig a free letter commuting with D (P e1 = e1, phase row 1 zero) of ODD order when available,
and also y0 sig y0^-1 free (y0 = D mu0), so that top 0 kernelises array 0 through BOTH merges:
V_0 = (c1 + c2 sig) q10 and F_00 = (c1 + c2 b') b (a two-term, generically invertible entry: the sigma-pair).
Array 1 is routed through merge 2 (x1 -> s1 -> u2 -> v -> t1).  Variant 'sig1': a11 random letter; 'sig0': a11 absent
(u_1 private to array 0).  Census of top-layer tags (Lemmas A', C), then the exact producer sweep (sweep2) and Singular.
Usage: sweep_sig.py variant ncores seed maxpool cap outdir"""
import sys, os, random, itertools, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sweep2 import *
import sweep2

def comm_letter(rng):
    """free letter commuting with D symplectically: P e1 = e1, B row/col 1 = 0"""
    while True:
        P = sp.gl4(rng)
        if P[0] != 1: 
            # force P e1 = e1: replace first column, keep invertibility
            P = [1] + P[1:]
            if len(sp.span(P)) != 16: continue
        Bm = [[0] * 4 for _ in range(4)]
        for i in range(1, 4):
            for j in range(i, 4):
                Bm[i][j] = Bm[j][i] = rng.randrange(2)
        S = sp.symp_from(P, Bm)
        return S, P, Bm

def order(M, maxo=64):
    A = M.copy()
    for k in range(1, maxo):
        if np.allclose(A, np.eye(16)): return k
        A = G.snap(A @ M)
    return None

def core_sig(rng, variant, tries=300):
    Dinv = np.conj(G.D1).T
    for _ in range(tries):
        R0 = route(rng)
        y0S = sp.comp(sp.D, R0['mu'][0])
        y0 = G.mul(G.D1, R0['mu'][1])
        found = None
        for _ in range(3000):
            S, P, Bm = comm_letter(rng)
            # y0 sig y0^-1 must stabilise L0
            C = sp.comp(y0S, sp.comp(S, sp.inv(y0S)))
            if sp.lag(C) != sp.L0: continue
            M = sp.to_matrix(P, Bm, rng.randrange(16), 0)
            g2 = G.snap(Dinv @ M @ G.D1 @ R0['g'][1])
            if not G.is_monomial(g2): continue
            if not G.is_monomial(G.snap(y0 @ M @ np.conj(y0).T)): continue
            o = order(M)
            if o is None or o == 1: continue
            found = (S, M, g2, o)
            if o % 2 == 1: break
        if found is None: continue
        S, M, g2, o = found
        g2S = sp.comp(sp.inv(sp.D), sp.comp(S, sp.comp(sp.D, R0['g'][0])))
        # array 1 through merge 2: b1, a21 with q21 = D g2 D a21 D b1 layer 3 and K-adjacent, then mu1
        for _ in range(400):
            b1 = rand_full(rng); s1 = sp.mv(b1[0], sp.I8)
            a21 = rand_full(rng); p21 = sp.mv(a21[0], s1)
            if sp.lam(p21) != 2: continue
            q21 = sp.mv(g2S, p21)
            if sp.lam(q21) == 3 and sp.dist(sp.lag(q21), sp.LK) == 1:
                mu1 = None
                for _ in range(3000):
                    m = rand_full(rng)
                    if sp.lag(sp.mv(m[0], q21)) == sp.LK: mu1 = m; break
                if mu1 is None: continue
                if sp.lag(sp.comp(sp.D, mu1[0])) == sp.lag(y0S): continue     # both tops in one frame: F singular
                R1 = dict(b=b1, a=a21, g=(g2S, g2), mu=mu1, s=s1, p=p21, q=q21)
                a20 = R0['a']
                for _ in range(200 if variant == 'sig1' else 1):
                    a11 = rand_full(rng) if variant == 'sig1' else None
                    tg = tags_ok(R0, R1, a11, a20)
                    # top 0 serves (1,0) and (2,0): the tag of (2,0) at top 0 must be 'K'
                    i20 = 1 if a11 is not None else 0
                    if tg[i20] != 'K': break
                    if not any(t in ('4', '3') for t in tg):
                        return R0, R1, a11, a20, tg, o
                break
    return None

if __name__ == '__main__':
    variant, ncores, seed, maxpool, cap, outdir = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed); tally = {}
    for c in range(ncores):
        t0 = time.time()
        cr = core_sig(rng, variant)
        if cr is None: print("core failed", flush=True); continue
        R0, R1, a11, a20, tg, o = cr
        mv0 = moves_of(R0, R1, a11, a20)
        sweep2.EXIST.clear()
        for mi, summ in enumerate(mv0):
            sweep2.EXIST[mi] = [(src, Mx) for (src, Mx) in summ]
        P0 = build(2, mv0)
        cells = []; sz = []
        for l0 in range(2):
            C = Closure(P0, l0, cap)
            if not C.ok: sz.append('cap'); continue
            sz.append((C.nr, C.nc))
            inv = {v: k for k, v in C.rows.items()}
            cells += [inv[r] for r in range(C.nr)]
        pool = producers(P0, R0, R1, cells)
        print(f"core {c}: sigma order {o}, tags {tg}, bare closures {sz}, exact producers {len(pool)}", flush=True)
        pool = pool[:maxpool]
        def expand(items):
            ex = []
            for it in items:
                if it[0] == 'pair': ex += [it[1], it[2]]
                else: ex.append(it)
            return ex
        subsets = [()] + [(i,) for i in range(len(pool))] + list(itertools.combinations(range(len(pool)), 2))
        res = {}
        for sub in subsets:
            mvs = moves_of(R0, R1, a11, a20, expand([pool[i] for i in sub]))
            P = build(2, mvs)
            scr, sizes = singular_script(P, cap=cap)
            if scr is None: v = 'CAP'
            else: v, out = decide(scr, os.path.join(outdir, f'{variant}_{seed}_{c}.sing'), timeout=90)
            res[v] = res.get(v, 0) + 1; tally[v] = tally.get(v, 0) + 1
            if v not in ('INCONSISTENT', 'CAP'):
                print(f"   subset {sub} sizes {sizes}: {v}", flush=True)
                if v == 'CONSISTENT':
                    pickle.dump({'moves': mvs}, open(os.path.join(outdir, f'CONS_{variant}_{seed}_{c}_{"_".join(map(str, sub))}.pkl'), 'wb'))
        print(f"core {c} ({variant}): pool {len(pool)}, {len(subsets)} skeletons: {res} ({time.time()-t0:.0f}s)", flush=True)
    print("TOTAL", tally, flush=True)
