"""topdiag.py (mw2-b5beyond, resumed 14 Oct ~14:45): WHY does every block-5 two-array skeleton stop at res_5 = 1.00?
For each skeleton: LM on the TOP layer only (d = m), then report per-array residuals and the rows that keep a residual
(coordinate l, layer, dK = distance of the term's frame from K's frame).  If the top system has rank 1 (the arrays only
reach K in a fixed combination), each array keeps 0.5.  Kinds:
  pure   : t_lemmaD's butterfly (single-letter starts and merge, pure multi-term trunk, ends with side terms)
  inject : the same with an input x_2 injected into the merge t1
  sepT   : v-ms-fano's separation inside the trunk: t1 = D(s1 + s2), t1' = D(s1 - sigma s2) (sigma = order-3 address
           map), both feed t2 (2 in 8 at block 5: s1 s2 t1 t1' t2 t3 e1 e2)
  plain  : two separate 5-chains (2 in 10) -- positive control, must reach res_5 = 0
usage: topdiag.py kind nskel seed"""
import sys, random, time
from b5eng import *
from lsearch2 import rand_free, dstep


def geo(E, rng, g, tries=3000):
    for _ in range(tries):
        b = rand_free(E, rng)
        if dstep(E, b, g) == 1:
            return b
    return None


def leads(E, P, gi):
    L = {}
    for (l, g) in P.gens[gi]:
        if l not in L or E.dK(g) < E.dK(L[l]):
            L[l] = g
    return L


def sigma3(E):
    m = E.m
    P = [1 << t for t in range(m)]
    P[0], P[1], P[2] = 1 << 1, 1 << 2, 1 << 0          # e0 -> e1 -> e2 -> e0
    return Letters(E).addr(P)


def build(kind, E, rng):
    P = Prog(E, 2); m = E.m
    if kind == 'plain':
        for l in range(2):
            prev = l
            for j in range(m):
                g = P.gens[prev][0][1] if prev >= 2 else E.ID
                b = geo(E, rng, g)
                prev = P.add_move([(prev, b)])
        return P
    b1 = geo(E, rng, E.ID); b2 = geo(E, rng, E.ID)
    s1 = P.add_move([(0, b1)]); s2 = P.add_move([(1, b2)])
    g1 = P.gens[s1][0][1]; g2 = P.gens[s2][0][1]
    a1 = geo(E, rng, g1); a2 = geo(E, rng, g2)
    if kind == 'sepT':
        sg = sigma3(E)
        t1 = P.add_move([(s1, a1), (s2, a2)])
        t1b = P.add_move([(s1, a1), (s2, -sg @ a2)])
        L = leads(E, P, t1); Lb = leads(E, P, t1b)
        bs = [(t1, geo(E, rng, L[0])), (t1, geo(E, rng, L[1])), (t1b, geo(E, rng, Lb[0])), (t1b, geo(E, rng, Lb[1]))]
        prev = P.add_move(bs)
        nt = 1
    else:
        summ = [(s1, a1), (s2, a2)]
        if kind == 'inject':
            summ.append((1, rand_free(E, rng)))
        prev = P.add_move(summ); t1 = prev
        nt = 2
    for k in range(nt):
        L = leads(E, P, prev)
        bs = [geo(E, rng, L[0]), geo(E, rng, L[1]), rand_free(E, rng)]
        prev = P.add_move([(prev, b) for b in bs])
    t3 = prev
    for j in range(2):
        L = leads(E, P, t3)
        summ = [(t3, geo(E, rng, L[0])), (t3, geo(E, rng, L[1]))]
        for src in [s1, s2, t1, t3 - 1, 0, 1]:
            summ.append((src, rand_free(E, rng)))
        P.add_move(summ)
    return P


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
            out.append((1.0, 0, 0)); continue
        beta = np.linalg.lstsq(A, t, rcond=1e-10)[0]
        r = A @ beta - t
        sv = np.linalg.svd(A, compute_uv=False)
        out.append((float(np.sum(np.abs(r) ** 2)), int(np.sum(sv > 1e-8 * sv[0])), nc))
    return out


def main():
    kind, nskel, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    rng = random.Random(seed); t0 = time.time()
    for s in range(nskel):
        E = Eng(5, seed=seed * 100 + s)
        try:
            P = build(kind, E, rng)
        except (KeyError, TypeError) as ex:
            print("skel %d: build failed %r" % (s, ex), flush=True); continue
        if not P.setup_final(cap=3000, dmin=5):
            print("skel %d: refused %s" % (s, P.why), flush=True); continue
        fi = P.free_idx(); nf = len(fi)
        best = (np.inf, None)
        for st in range(4):
            x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
            x, c = P.lm(5, x0, fi, iters=80)
            if c < best[0]:
                best = (c, x)
        pa = per_array(P, P.full_alpha(best[1]), 5)
        print("skel %d %s: terms %s closure %s | res_5=%.3e per array (res, rank, cols) %s (%.0fs)" % (
            s, kind, [len(g) for g in P.gens][2:], [(f[0], f[1]) for f in P.fin], best[0],
            [("%.3f" % a, b, c) for (a, b, c) in pa], time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
