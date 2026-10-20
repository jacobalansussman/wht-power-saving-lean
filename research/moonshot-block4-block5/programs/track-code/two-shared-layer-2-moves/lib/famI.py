"""famI.py -- w3-asym: closure-guided ("frame-matched") injections on the reduced BRIEF3 shape, judged on the FULL
exact-structure residual (closure test, all f at once).
Core (single letters): starts s_l = D b_l x_l; merge u = D(sum a_l s_l) with the clean split; end l: v_l = D(ell_l u),
t_l = D(m_l v_l) with D m_l D ell_l p_l in the K frame.  NOTE.md's reduction: output l's readout uses the top moves only
through t'_l (here t_l), so its junk w_l p_j x_j (j != l) and everything it drags in must be cancelled by lower generators.
Greedy loop: for every output, compute the closure and the least-squares residual; take the rows (cells) with the largest
residual; for each such cell look for injections (target move in {u, v_e, t_e}, source an earlier generator, engine letter b
sampled) that put a term of the source into that cell's frame (directly, or after the end's next letter); add the best one
(min residual after LM over all coefficients); repeat while the move set allows.  Records every step.
Usage: famI.py W nskel seed maxinj nsamp outprefix"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from famP import *

def lamg(g):
    return lam(G.MATS[g])

def core(rng, W, tries=300):
    for _ in range(tries):
        bs = [prog_letter(rng, [G.ID]) for l in range(W)]
        hs = [G.intern(G.mul(G.D1, b)) for b in bs]
        as_ = [prog_letter(rng, [h]) for h in hs]
        ps = [G.mul(G.D1, a, G.MATS[h]) for a, h in zip(as_, hs)]
        if all(lam(G.mul(G.K, np.conj(ps[l]).T, ps[j])) <= 3 for l in range(W) for j in range(W) if l != j):
            break
    else:
        return None
    pids = [G.intern(p) for p in ps]
    rec = [[(l, bs[l])] for l in range(W)]                    # moves 0..W-1: starts (generator index W + l)
    rec.append([(W + l, as_[l]) for l in range(W)])             # merge u: generator 2W
    u = 2 * W
    for l in range(W):
        path = midpaths(rng, pids[l], nsamp=400)
        if not path: return None
        ell, m = path[0]
        rec.append([(u, ell)])                                  # v_l: generator 2W+1+2l
        rec.append([(u + 1 + 2 * l, m)])                        # t_l: generator 2W+2+2l
    return rec

def build_rec(W, rec):
    P = Prog(W)
    for summ in rec:
        P.add_move(summ)
    return P

def full_res(P, starts, nrng, cap=4000, x0=None):
    """LM over all parameters of the summed closure residual of all outputs; returns (res, par, per-output, closures)"""
    Cs = [Closure(P, l0, cap) for l0 in range(P.W)]
    if not all(C.ok for C in Cs):
        return None
    k = P.np
    def fun(x):
        par = x[:k] + 1j * x[k:]
        v = np.concatenate([C.resid_vec(par) for C in Cs])
        return np.concatenate([v.real, v.imag])
    best = (np.inf, None)
    if k == 0:
        par = np.zeros(0, dtype=complex)
        r = [C.residual(par) for C in Cs]
        return sum(r), par, r, Cs
    for st in range(starts):
        xs = nrng.standard_normal(2 * k) if (x0 is None or st > 0) else x0
        r = least_squares(fun, xs, method='trf', max_nfev=60 * (2 * k + 1))
        if 2 * r.cost < best[0]: best = (2 * r.cost, r.x)
    par = best[1][:k] + 1j * best[1][k:]
    return best[0], par, [C.residual(par) for C in Cs], Cs

def bad_cells(C, par, top=4):
    A = C.matrix(par); t = np.zeros(C.nr, dtype=complex); t[0] = 1
    beta = np.linalg.lstsq(A, t, rcond=None)[0]; r = np.abs(A @ beta - t) ** 2
    inv = {v: k for k, v in C.rows.items()}
    order = np.argsort(-r)
    return [inv[i] for i in order[:top] if r[i] > 1e-10]

def candidates(P, rec, W, cell, rng, nsamp):
    """injections putting a term into the frame of cell (l, g): (move index in rec, src, b, how)"""
    l, g = cell
    fk = G.frm(g)
    out = []
    nmoves = len(rec)
    for mi in range(W, nmoves):                 # target move (not a start)
        gi = W + mi                             # its generator index
        # the end letter applied after this move (if mi is a v-move, its t follows with letter m)
        nxt = None
        for mj in range(mi + 1, nmoves):
            if len(rec[mj]) >= 1 and rec[mj][0][0] == gi:
                nxt = rec[mj][0][1]; break
        for src in range(gi):                   # any earlier generator (inputs 0..W-1, starts W.., ...)
            hs = [h for (ll, h) in P.terms[src] if ll == l]
            if not hs: continue
            for _ in range(nsamp):
                b = rand_engine(rng)
                for h in hs:
                    t1 = G.intern(G.mul(G.D1, b, G.MATS[h]))
                    if G.frm(t1) == fk:
                        out.append((mi, src, b, 'direct')); break
                    if nxt is not None:
                        t2 = G.intern(G.mul(G.D1, nxt, G.MATS[t1]))
                        if G.frm(t2) == fk:
                            out.append((mi, src, b, 'next')); break
    # dedupe by (mi, src, how) keeping a few
    seen = {}; ded = []
    for c in out:
        key = (c[0], c[1], c[3])
        if seen.get(key, 0) < 2:
            seen[key] = seen.get(key, 0) + 1; ded.append(c)
    return ded

def run(W, nsk, seed, maxinj, nsamp, outp, corefn):
    rng = random.Random(seed); nrng = np.random.default_rng(seed)
    for s in range(nsk):
        t0 = time.time()
        rec = corefn(rng, W)
        if rec is None:
            print(f"skel {s}: core failed", flush=True); continue
        P = build_rec(W, rec)
        fr = full_res(P, 1, nrng)
        print(f"skel {s}: core residual {fr[0]:.4f} per output {np.round(fr[2], 4)} closures {[(C.nr, C.nc) for C in fr[3]]}", flush=True)
        cur = (fr[0], rec, fr)
        for step in range(maxinj):
            Pc = build_rec(W, cur[1]); res0, par0, per0, Cs = cur[2]
            cells = []
            for C in Cs:
                cells += bad_cells(C, par0, top=3)
            cands = []
            for cell in cells:
                cands += candidates(Pc, cur[1], W, cell, rng, nsamp)
            if not cands:
                print(f"   step {step}: no candidate injections for cells {[(c[0], lamg(c[1])) for c in cells]}", flush=True); break
            best = None
            for (mi, src, b, how) in cands[:24]:
                rec2 = [list(m) for m in cur[1]]; rec2[mi] = rec2[mi] + [(src, b)]
                P2 = build_rec(W, rec2)
                fr2 = full_res(P2, 2, nrng)
                if fr2 is None: continue
                if best is None or fr2[0] < best[0]:
                    best = (fr2[0], rec2, fr2, (mi, src, how))
            if best is None or best[0] > cur[0] - 1e-6:
                print(f"   step {step}: no improving injection among {len(cands)} (best {best[0] if best else None})", flush=True); break
            cur = best[:3]
            print(f"   step {step}: + injection move {best[3][0]} src g{best[3][1]} ({best[3][2]}): residual {best[0]:.5f} per output {np.round(best[2][2], 4)} closures {[(C.nr, C.nc) for C in best[2][3]]}", flush=True)
            if best[0] < 1e-9:
                pickle.dump({'W': W, 'rec': cur[1], 'par': best[2][1]}, open(f"{outp}_hit{s}.pkl", 'wb'))
                print("   HIT (numerical): saved", flush=True); break
        print(f"skel {s} done: final residual {cur[0]:.5f} ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6], core)
