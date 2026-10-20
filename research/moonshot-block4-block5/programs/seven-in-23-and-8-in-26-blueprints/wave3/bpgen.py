"""bpgen.py -- w3-blueprint: random STRUCTURED skeletons for the blueprint orders (moonshot wave 3, 14 Oct 2026).

shape: string over {F, O}: F = tau-fixed move (orbit-sum content, a tau-fixed direction), O = free orbit (W moves).
Contents: up to `menu` summands (source generator, free letter); sources = raw inputs (own array first), earlier
generators (any orbit member); letters = progress letters of a lead term of the source (the deepest terms first) or an
ALIGNED letter putting a term of another source onto the first summand's new term (a collision).  The first free orbit
always takes its own raw input (mirror-D2 needs W independent co-top routes).
Score: frame-sum relaxation residual for array 0 (LM over the content coefficients).  By the character splitting
(NOTE.md) residual = sum over characters chi of H of the chi-share (each <= 1/W); k/W means k characters fail.
Skeletons with relaxation < 1e-8 get the FULL exact-structure test (closure for array 0 + VP-LM).
Usage: bpgen.py W H shape menu nskel seed0 cap"""
import sys, time
from symsolve import *
from bp import frame_relax, prog
import os
VERB = bool(os.environ.get('BPV'))
TWIST = bool(os.environ.get('BPTW'))

def build(W, H, shape, menu, seed, R=None):
    sk = Skel(4, W, H, seed, R)
    G, rng, R = sk.G, sk.rng, sk.R
    fd = fixed_dirs(R, rng)
    first_free = True
    for pos, tok in enumerate(shape):
        n = len(sk.gens)
        # candidate (source, term) pairs: terms of array 0 in any generator, deepest (smallest distance) first
        pool = []
        for i in range(n):
            for (l, g) in sk.gens[i]:
                dd = G.dist(g)
                if dd >= 1:
                    pool.append((dd, i, l, g))
        rng.shuffle(pool)
        pool.sort(key=lambda t: (t[0] + (0 if t[2] == 0 else 0.5) + rng.random() * 1.5))
        summ = []
        if tok == 'F' and TWIST and pos < len(shape) - 1:
            # SELECTIVE shared move (NOTE.md 5b): exact orbit sum, x-direction, content letter b with (b g)^-1 c (b g)
            # outside K^-1 B K (twisting), and progress dist(d b g) = dist(g) - 1
            c = R.Hc[1][0]; Kc = np.conj(G.K).T
            found = None
            for (dd, i, l, g) in pool[:20]:
                for d in fd:
                    for _ in range(300):
                        b = rand_free_b(G, rng, 4)
                        bg = G.mul(b, G.MATS[g])
                        if G.dist(G.intern(G.mul(d, bg))) != dd - 1:
                            continue
                        if G.is_monomial(G.mul(G.K, np.conj(bg).T, c, bg, Kc)):
                            continue
                        found = (d, i, b, g); break
                    if found: break
                if found: break
            if found is None:
                print(f"  build: no twisting progress letter at token {pos}", flush=True); return None
            fkind = 'S'
            Dm, i, b, g = found
            summ.append((i, b))
        elif tok == 'F':
            found = None
            for (dd, i, l, g) in pool[:60]:
                for d in rng.sample(fd, len(fd)):
                    pl = prog(G, g, rng, Dmat=d, want=8)
                    if pl:
                        found = (d, i, rng.choice(pl), g); break
                if found:
                    break
            fkind = 'S'
            if found is None or rng.random() < 0.5:
                # c-twisted orbit sum: absorbed move in any direction
                for (dd, i, l, g) in pool[:60]:
                    pl = prog(G, g, rng, want=4)
                    if pl:
                        found = (G.D1, i, rng.choice(pl), g); fkind = 'T'; break
            if found is None:
                print(f"  build: no fixed move progresses at token {pos} ({shape})", flush=True); return None
            Dm, i, b, g = found
            summ.append((i, b))
        else:
            Dm = G.D1
            if first_free:
                cands = [(4, 0, 0, G.ID)] + [t for t in pool if t[1] != 0]
                first_free = False
            else:
                cands = pool
            found = None
            for (dd, i, l, g) in cands[:12]:
                pl = prog(G, g, rng, want=4)
                if pl:
                    found = (i, rng.choice(pl), g); break
            if found is None:
                print(f"  build: no progress letter at token {pos} ({shape}); pool dists {sorted(set(t[0] for t in cands))}", flush=True); return None
            i, b, g = found
            summ.append((i, b))
        T = G.intern(G.mul(Dm, summ[0][1], G.MATS[found[-1]]))
        # more summands: half aligned (collision onto T), half progress letters of other lead terms
        tries = 0
        while len(summ) < menu and tries < 40:
            tries += 1
            if rng.random() < 0.5:
                (dd, i, l, g) = rng.choice(pool) if pool else (4, 0, 0, G.ID)
                if l != 0 and rng.random() < 0.5:
                    continue
                b = G.mul(np.conj(Dm).T, G.MATS[T], np.conj(G.MATS[g]).T)
                if G.is_monomial(b) and (i, G.key(b)) not in [(s, G.key(bb)) for (s, bb) in summ]:
                    summ.append((i, b))
            else:
                src = rng.choice(pool + [(4, 0, 0, G.ID)] * 3)
                (dd, i, l, g) = src
                pl = prog(G, g, rng, Dmat=Dm, want=3)
                if pl:
                    summ.append((i, rng.choice(pl)))
        t1 = time.time()
        sk.add_move(fkind if tok == 'F' else 'F', summ, Dm if tok == 'F' else None)
        if VERB:
            print(f"   token {pos} {tok}: {len(summ)} summands, add_move {time.time()-t1:.1f}s, rep terms {sk.moves[-1]['n']}", flush=True)
    cost = sum(1 + len(mv["images"]) for mv in sk.moves)
    return sk, cost

def main():
    W, H, shape, menu, nskel, seed0, cap = int(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
    tlimit = float(sys.argv[8]) if len(sys.argv) > 8 else 1e9
    t0 = time.time(); res = []
    R = Run(4, W, H, seed0)
    for s in range(nskel):
        if time.time() - t0 > tlimit:
            break
        seed = seed0 + s
        try:
            R.rng = random.Random(seed)
            r = build(W, H, shape, menu, seed, R)
        except Exception as e:
            print(f"seed {seed}: build error {type(e).__name__} {e}", flush=True); continue
        if r is None:
            print(f"seed {seed}: build failed", flush=True); continue
        sk, cost = r
        orbs = [1 + len(mv['images']) for mv in sk.moves]
        tb = time.time()
        try:
            rel = frame_relax(sk, starts=2)
        except Exception as e:
            print(f"seed {seed}: relax error {e}", flush=True); continue
        k = None if rel is None else rel * W
        if VERB:
            print(f"   relax took {time.time()-tb:.1f}s", flush=True)
        print(f"seed {seed}: cost {cost} orbits {orbs} terms {[mv['n'] for mv in sk.moves]} relax res^2 {rel} (= {k} chars) ({time.time()-t0:.0f}s)", flush=True)
        res.append((rel if rel is not None else 9.0, seed, cost))
        if rel is not None and rel < 1e-8 and cost == (23 if W == 7 else 26):
            Wsave = sk.W; sk.W = 1
            ok = sk.setup_final(cap); sk.W = Wsave
            if not ok:
                print(f"seed {seed}: FULL closure cap {cap}", flush=True); continue
            prep_fast(sk)
            N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
            best = np.inf
            for st in range(2):
                x0 = np.array([sk.rng.gauss(0, 1) for _ in range(2 * nf)])
                x, c = vp_lm(sk, N, free_idx, x0)
                best = min(best, c)
                if c < 1e-14:
                    np.save(f"bpgen_sol_{W}_{H}_{shape}_{seed}.npy", x)
            print(f"seed {seed}: FULL dims {[(nr, nc) for (nr, nc, E) in sk.final]} unknowns {nf} best res^2 {best:.4e}", flush=True)
    res.sort()
    print("summary best (relax, seed, cost):", res[:10], flush=True)
    from collections import Counter
    print("relax histogram (chars failing):", sorted(Counter(round(r * W, 2) for (r, s, c) in res).items()), flush=True)

if __name__ == "__main__":
    main()
