"""sweep.py -- w5-twoshared: exact sweep of the 2+2+1+2 shape (W = 2, block 4, engine letters, every coefficient free).
Generators: 0 x0, 1 x1, 2 s0 = D b0 x0, 3 s1 = D b1 x1, 4 u1 = D(a10 s0 + a11 s1), 5 u2 = D(a20 s0 + a21 s1),
6 v = D(g1 u1 + g2 u2), 7 t0 = D(mu0 v + ...), 8 t1 = D(mu1 v + ...).  Routes first (F diagonal): x0->s0->u1->v->t0,
x1->s1->u2->v->t1.  variant 'priv': a11 = a20 = 0 (layer-2 private, v the only shared move); 'cross': both merges carry
both arrays, cross letters re-sampled until the top-layer census conditions hold (no non-kernel layer-4 junk, no
layer-3 junk outside V's frames: Lemmas A', C).  Injection pool = closure-guided landing injections: for an uncancelled
closure cell (array j, frame Phi) of the bare skeleton, letters h with frame(D h src_term) = Phi for sources
s_j, u_1, u_2 into t_0 / t_1 / v, and x_j into u_1 / u_2 / v (found in the symplectic picture, lifted with random
Pauli parts).  Every subset of size <= 2 of the pool (plus the bare skeleton) is decided by Singular.
Usage: sweep.py variant ncores seed npool cap outdir"""
import sys, os, random, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sp
from exact import *

def lift(rng, L):
    S, P, Bm = L
    return S, sp.to_matrix(P, Bm, rng.randrange(16), rng.randrange(16))

def rand_full(rng):
    return lift(rng, sp.letter_full(rng))

def route(rng):
    """route letters (symplectic, exact): b, a, g, mu with D mu D g D a D b in the kernel frame"""
    while True:
        b = rand_full(rng); s = sp.mv(b[0], sp.I8)
        for _ in range(50):
            a = rand_full(rng); p = sp.mv(a[0], s)
            if sp.lam(p) == 2: break
        else: continue
        for _ in range(600):
            g = rand_full(rng); q = sp.mv(g[0], p)
            if sp.lam(q) == 3 and sp.dist(sp.lag(q), sp.LK) == 1:
                for _ in range(3000):
                    mu = rand_full(rng)
                    if sp.lag(sp.mv(mu[0], q)) == sp.LK:
                        return dict(b=b, a=a, g=g, mu=mu, s=s, p=p, q=q)
                break

def tags_ok(R0, R1, a11, a20):
    Q = {(1, 0): R0['q'], (2, 1): R1['q']}
    if a11 is not None:
        Q[(1, 1)] = sp.mv(R0['g'][0], sp.mv(a11[0], R1['s']))
    if a20 is not None:
        Q[(2, 0)] = sp.mv(R1['g'][0], sp.mv(a20[0], R0['s']))
    Y = {0: sp.comp(sp.D, R0['mu'][0]), 1: sp.comp(sp.D, R1['mu'][0])}
    tags = []
    for e, rt in ((0, (1, 0)), (1, (2, 1))):
        for kj in sorted(Q):
            if kj == rt: continue
            J = sp.comp(Y[e], Q[kj]); L = sp.lam(J); LJ = sp.lag(J)
            if L == 4: tags.append('K' if LJ == sp.LK else '4')
            elif L == 3:
                tags.append('3V' if any(sp.lag(Q[(k2, j2)]) == LJ for (k2, j2) in Q if j2 == kj[1]) else '3')
            else: tags.append(str(L))
    return tags

def core(rng, variant, tries=400):
    for _ in range(tries):
        R0 = route(rng); R1 = route(rng)
        if variant == 'priv':
            a11 = a20 = None
            tg = tags_ok(R0, R1, None, None)
            if any(t in ('4', '3', 'K') for t in tg): continue
            return R0, R1, None, None, tg
        for _ in range(200):
            a11 = rand_full(rng); a20 = rand_full(rng)
            tg = tags_ok(R0, R1, a11, a20)
            if not any(t in ('4', '3', 'K') for t in tg):
                return R0, R1, a11, a20, tg
    return None

def moves_of(R0, R1, a11, a20, extra=()):
    u1 = [(2, R0['a'][1])] + ([(3, a11[1])] if a11 is not None else [])
    u2 = [(3, R1['a'][1])] + ([(2, a20[1])] if a20 is not None else [])
    mv = [[(0, R0['b'][1])], [(1, R1['b'][1])], u1, u2, [(4, R0['g'][1]), (5, R1['g'][1])],
          [(6, R0['mu'][1])], [(6, R1['mu'][1])]]
    for (mi, src, M) in extra:
        mv[mi] = mv[mi] + [(src, M)]
    return mv

def symp_terms(P, i):
    """symplectic images of generator i's terms: list of (l, S) (recomputed from exact matrices via the term words)"""
    return None

def landing(rng, P, cell, nsamp=60):
    """injections (move index, src, matrix) putting a term of src exactly into the frame of cell (l, gid).
    Moves: 2 u1, 3 u2, 4 v, 5 t0, 6 t1 (indices into the move list; generator = move + 2)."""
    l, gid = cell
    fk = G.frm(gid)
    out = []
    for mi in (2, 3, 4, 5, 6):
        gi = mi + 2
        for src in range(gi):
            hs = [h for (ll, h) in P.terms[src] if ll == l]
            for h in hs:
                for _ in range(nsamp):
                    S, Pp, Bm = sp.letter_full(rng)
                    M = sp.to_matrix(Pp, Bm, rng.randrange(16), rng.randrange(16))
                    t1 = G.mul(G.D1, M, G.MATS[h])
                    if G.frame_key(t1) == fk:
                        out.append((mi, src, M)); break
    return out

if __name__ == '__main__':
    variant, ncores, seed, npool, cap, outdir = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed); tally = {}
    for c in range(ncores):
        t0 = time.time()
        cr = core(rng, variant)
        if cr is None: print("core failed", flush=True); continue
        R0, R1, a11, a20, tg = cr
        P0 = build(2, moves_of(R0, R1, a11, a20))
        cells = []
        for l0 in range(2):
            C = Closure(P0, l0, cap)
            if not C.ok: continue
            inv = {v: k for k, v in C.rows.items()}
            # rows of the closure other than the target: cells that must be cancelled
            cells += [inv[r] for r in range(1, C.nr)]
        print(f"core {c}: tags {tg}, bare closures {[Closure(P0, l0, cap).nr if Closure(P0, l0, cap).ok else 'cap' for l0 in range(2)]}, {len(cells)} cells", flush=True)
        rng.shuffle(cells)
        pool = []
        for cell in cells[:8]:
            pool += landing(rng, P0, cell, nsamp=40)
            if len(pool) >= npool: break
        pool = pool[:npool]
        subsets = [()] + [(i,) for i in range(len(pool))] + list(itertools.combinations(range(len(pool)), 2))
        res = {}
        for sub in subsets:
            P = build(2, moves_of(R0, R1, a11, a20, [pool[i] for i in sub]))
            scr, sizes = singular_script(P, cap=cap)
            if scr is None: v = 'CAP'
            else: v, out = decide(scr, os.path.join(outdir, f'{variant}_{seed}_{c}.sing'), timeout=60)
            res[v] = res.get(v, 0) + 1; tally[v] = tally.get(v, 0) + 1
            if v not in ('INCONSISTENT', 'CAP'):
                print(f"   subset {[(pool[i][0], pool[i][1]) for i in sub]} sizes {sizes}: {v}", flush=True)
                if v == 'CONSISTENT':
                    import pickle
                    pickle.dump({'moves': moves_of(R0, R1, a11, a20, [pool[i] for i in sub])}, open(os.path.join(outdir, f'CONS_{variant}_{seed}_{c}_{len(sub)}.pkl'), 'wb'))
        print(f"core {c} ({variant}): pool {len(pool)}, {len(subsets)} skeletons: {res} ({time.time()-t0:.0f}s)", flush=True)
    print("TOTAL", tally, flush=True)
