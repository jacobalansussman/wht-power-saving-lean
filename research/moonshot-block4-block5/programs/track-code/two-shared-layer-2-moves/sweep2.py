"""sweep2.py -- w5-twoshared: as sweep.py, but the injection pool is built by the EXACT PRODUCER RULE: for every closure
cell of the bare skeleton that holds a term J of a top t_e (array j), every single injection whose new term inside t_e
equals J up to a scalar is computed directly (no sampling):
   src -> t_e   with letter M = D^-1 J h^-1,
   src -> v     with letter M = (D mu_e D)^-1 J h^-1        (src in {x, s, u_1, u_2}),
   src -> u_k   with letter M = (D mu_e D g_k D)^-1 J h^-1  (src in {x, s}),
for every term h of src on array j; kept iff M is a monomial with entries in {1, i, -1, -i} (an engine/widened letter).
Then the bare skeleton, every single producer and every pair are decided by Singular (every coefficient free).
Usage: sweep2.py variant ncores seed maxpool cap outdir"""
import sys, os, random, itertools, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sweep import *

def is_Z4_mono(M):
    return G.is_monomial(M)

def lam_ex(M):
    return int(round(np.log2((np.abs(M) > 1e-9).sum(axis=1).max())))

EXIST = {}
def existing(mi, src, M):
    for (s2, M2) in EXIST.get(mi, []):
        if s2 == src:
            A = M2; i = np.unravel_index(np.argmax(np.abs(A)), A.shape); z = M[i] / A[i]
            if abs(abs(z) - 1) < 1e-9 and np.allclose(M, z * A, atol=1e-9): return True
    return False

def factor1(X, tries=4000, rng=random.Random(99)):
    # X at layer 1: find free h with h^-1 X in D B' (then kappa = D^-1 h^-1 X)
    Dinv = np.conj(G.D1).T
    for _ in range(tries):
        S, Pp, Bm = sp.letter_full(rng)
        h = sp.to_matrix(Pp, Bm, rng.randrange(16), rng.randrange(16))
        kap = G.snap(Dinv @ np.conj(h).T @ X)
        if G.is_monomial(kap): return h, kap
    return None

def producers(P, R0, R1, cells):
    Dinv = np.conj(G.D1).T
    mu = {5: R0['mu'][1], 6: R1['mu'][1]}          # move index -> mu letter of that top
    gk = {2: R0['g'][1], 3: R1['g'][1]}           # u_1 (move 2) via g1, u_2 (move 3) via g2
    out = []; seen = set()
    for (l, gid) in cells:
        fk = G.frm(gid)
        for tm in (5, 6):
            tg = tm + 2
            for (ll, h) in P.terms[tg]:
                if ll != l or G.frm(h) != fk: continue
                J = G.MATS[h]
                Dm = G.mul(G.D1, mu[tm])
                targets = [(tm, J, range(tg - 2))]                                    # into t_e: sources before v... any earlier gen except v itself is fine
                targets.append((4, G.mul(np.conj(Dm).T, J), range(6)))               # into v: sources x, s, u
                for um in (2, 3):
                    pre = G.mul(Dm, G.D1, gk[um])
                    targets.append((um, G.mul(np.conj(pre).T, J), range(4)))          # into u_k: sources x, s
                for (mi, Jt, srcs) in targets:
                    for src in srcs:
                        for (l2, h2) in P.terms[src]:
                            if l2 != l: continue
                            M = G.snap(Dinv @ Jt @ np.conj(G.MATS[h2]).T)
                            if is_Z4_mono(M):
                                key = (mi, src, G.key(M))
                                if key not in seen and not existing(mi, src, M):
                                    seen.add(key); out.append((mi, src, M))
                # (C4) pair: x_l -> u_k (kappa) and u_k -> t_e (h) with D h D kappa = J  <=>  D^-1 J = h D kappa (layer 1)
                X = G.mul(Dinv, J)
                if lam_ex(X) == 1:
                    for um in (2, 3):
                        hk = factor1(X)
                        if hk is not None:
                            h, kap = hk
                            out.append(('pair', (tm, um + 2, h), (um, l, kap)))
    return out

if __name__ == '__main__':
    variant, ncores, seed, maxpool, cap, outdir = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed); tally = {}; ptally = []
    for c in range(ncores):
        t0 = time.time()
        cr = core(rng, variant)
        if cr is None: print("core failed", flush=True); continue
        R0, R1, a11, a20, tg = cr
        mv0 = moves_of(R0, R1, a11, a20)
        EXIST.clear()
        for mi, summ in enumerate(mv0):
            EXIST[mi] = [(src, M) for (src, M) in summ]
        P0 = build(2, mv0)
        cells = []; sz = []
        for l0 in range(2):
            C = Closure(P0, l0, cap)
            if not C.ok: sz.append('cap'); continue
            sz.append((C.nr, C.nc))
            inv = {v: k for k, v in C.rows.items()}
            cells += [inv[r] for r in range(C.nr)]
        pool = producers(P0, R0, R1, cells)
        ptally.append(len(pool))
        print(f"core {c}: tags {tg}, bare closures {sz}, {len(cells)} cells, exact producers {len(pool)}: {[(p[0], p[1]) if p[0] != 'pair' else ('pair', p[1][0], p[2][0]) for p in pool]}", flush=True)
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
                print(f"   subset {[(pool[i][0], pool[i][1]) if pool[i][0] != 'pair' else 'pair' for i in sub]} sizes {sizes}: {v}", flush=True)
                if v == 'CONSISTENT':
                    pickle.dump({'moves': mvs}, open(os.path.join(outdir, f'CONS_{variant}_{seed}_{c}_{"_".join(map(str, sub))}.pkl'), 'wb'))
        print(f"core {c} ({variant}): pool {len(pool)}, {len(subsets)} skeletons: {res} ({time.time()-t0:.0f}s)", flush=True)
    print("TOTAL", tally, "producer counts", ptally, flush=True)
