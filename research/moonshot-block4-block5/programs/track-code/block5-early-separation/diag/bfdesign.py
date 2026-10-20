"""bfdesign.py (mw2-b5beyond, resumed 14 Oct ~14:55): FRAME-EXACT butterfly designs at block 5 (2 in 7 and 2 in 8),
to see at which layer the surviving butterfly family fails, and whether input injections into the early trunk (the
only escape that Lemma Delta leaves, PROOF.md section 3b) or v-ms-fano's order-3 separation change that.

Design (no top junk by construction):
  s_j = D b_j x_j                (geodesic first letters)
  t1  = D(a1 s1 + a2 s2 [+ inj])  a_j progress for s_j's lead
  [sep: t1' = D(a1 s1 - sigma a2 s2 [+ inj]), sigma = order-3 coordinate cycle: 2 in 8]
  t2  = D(c t1 [+ c' t1'] [+ inj])  c JOINT progress letter for both leads (single letter)
  t3  = D(c3 t2)                  joint progress
  e_k = D(split_k t3 + side terms)  split_k: progress for array k's lead only
inj: n_inj letters on EACH input (x1 and x2), random free letters or sigma-twisted; side: n_side random side terms at the
ends.  Layered LM: res_d for d = 5, 4, ..., 1 (stop at first failure); per-array residuals at the failing layer.
usage: bfdesign.py variant nskel seed n_inj n_side      variant in {pure, inj1, inj2, inj12, sep, sep_inj12}
  inj1: inject into t1 only; inj2: into t2 only; inj12: both."""
import sys, random, time
from b5eng import *
from lsearch2 import rand_free, dstep, main_letter, leads as leads2


def leads(E, P, gi):
    L = leads2(E, P, gi)
    return {l: L[l][1] for l in L}


def sigma3(E):
    m = E.m
    Pm = [1 << t for t in range(m)]
    Pm[0], Pm[1], Pm[2] = 1 << 1, 1 << 2, 1 << 0
    return Letters(E).addr(Pm)


def geo(E, rng, g, tries=3000):
    for _ in range(tries):
        b = rand_free(E, rng)
        if dstep(E, b, g) == 1:
            return b
    return None


def build(variant, E, rng, n_inj, n_side):
    P = Prog(E, 2); tags = []
    sg = sigma3(E)
    b1 = geo(E, rng, E.ID); b2 = geo(E, rng, E.ID)
    s1 = P.add_move([(0, b1)]); s2 = P.add_move([(1, b2)])
    a1 = geo(E, rng, P.gens[s1][0][1]); a2 = geo(E, rng, P.gens[s2][0][1])

    def inj():
        out = []
        for l in (0, 1):
            for _ in range(n_inj):
                b = rand_free(E, rng)
                if rng.random() < 0.5:
                    b = sg @ b
                out.append((l, b))
        return out
    i1 = 'inj1' in variant or 'inj12' in variant
    i2 = 'inj2' in variant or 'inj12' in variant
    t1 = P.add_move([(s1, a1), (s2, a2)] + (inj() if i1 else []))
    srcs = [t1]
    if variant.startswith('sep'):
        t1b = P.add_move([(s1, a1), (s2, -sg @ a2)] + (inj() if i1 else []))
        srcs.append(t1b)
    gl = leads(E, P, t1)
    c, tg = main_letter(E, rng, gl, None); tags.append(tg)
    summ = [(t1, c)]
    if variant.startswith('sep'):
        glb = leads(E, P, srcs[1])
        cb, tg = main_letter(E, rng, glb, None); tags.append(tg)
        summ.append((srcs[1], cb))
    t2 = P.add_move(summ + (inj() if i2 else []))
    gl = leads(E, P, t2)
    c3, tg = main_letter(E, rng, gl, None); tags.append(tg)
    t3 = P.add_move([(t3src, c3) for t3src in [t2]])
    gl = leads(E, P, t3)
    ends = []
    for k in (0, 1):
        sp, tg = main_letter(E, rng, gl, k); tags.append(tg)
        summ = [(t3, sp)]
        cand = [0, 1, s1, s2, t1, t2] + srcs[1:]
        for _ in range(n_side):
            summ.append((cand[rng.randrange(len(cand))], rand_free(E, rng)))
        ends.append(P.add_move(summ))
    return P, tags


def per_array(P, alpha, d):
    fl = P.flat(alpha); out = []
    for (nr, nc, R, C, F, lay) in P.fin:
        keep = lay >= d
        idx = -np.ones(nr, dtype=np.int64); idx[keep] = np.arange(keep.sum())
        sel = keep[R]
        R2, C2, F2 = idx[R[sel]], C[sel], F[sel]
        nr2 = int(keep.sum())
        A = np.zeros((nr2, nc), dtype=complex)
        np.add.at(A, (R2, C2), fl[F2])
        t = np.zeros(nr2, dtype=complex); t[idx[0]] = 1
        if nc == 0:
            out.append(1.0); continue
        beta = np.linalg.lstsq(A, t, rcond=1e-10)[0]
        out.append(float(np.sum(np.abs(A @ beta - t) ** 2)))
    return out


def main():
    variant, nskel, seed, n_inj, n_side = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    rng = random.Random(seed); t0 = time.time(); hist = {}
    for s in range(nskel):
        E = Eng(5, seed=seed * 100 + s)
        P, tags = build(variant, E, rng, n_inj, n_side)
        if not P.setup_final(cap=3000, maxcells=4e7, maxterms=4000, dmin=1):
            print("skel %d: refused %s" % (s, P.why), flush=True); continue
        fi = P.free_idx(); nf = len(fi); res = {}; pa = None
        for d in range(5, 0, -1):
            best = (np.inf, None)
            for st in range(3):
                x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
                x, cst = P.lm(d, x0, fi, iters=80)
                if cst < best[0]:
                    best = (cst, x)
                if best[0] < 1e-20:
                    break
            res[d] = best[0]
            if best[0] > 1e-10:
                pa = per_array(P, P.full_alpha(best[1]), d)
                break
        deep = min([d for d in res if res[d] < 1e-10], default=6)
        hist[deep] = hist.get(deep, 0) + 1
        print("skel %d %s: tags %s unknowns %d terms %s closure %s | %s | per-array at fail %s (%.0fs)" % (
            s, variant, tags, nf, [len(g) for g in P.gens][2:], [(f[0], f[1]) for f in P.fin],
            " ".join("d%d=%.1e" % (d, res[d]) for d in sorted(res, reverse=True)),
            None if pa is None else ["%.3f" % v for v in pa], time.time() - t0), flush=True)
    print("SUMMARY %s n_inj %d n_side %d: deepest layer passed (6 = none): %s" % (variant, n_inj, n_side, sorted(hist.items())), flush=True)


if __name__ == "__main__":
    main()
