"""lsearch.py (mw2-b5beyond): layered numerical search over skeletons of a given dependency shape.

usage: lsearch.py SHAPE W m menu s_e s_in starts seed nskel [cap] [dlist]
  SHAPE: name in SHAPES or 'chainT' (a T-move chain); covers (i, j) = output of move i is in the content of move j.
  menu: 'singer' (eta^k x X-Pauli), 'levi' (odd-order Levi elements, powers, x Pauli x Clifford phase),
        'rand' (random affine x Clifford phase x Pauli; no structure; control), 'levi1p' (1 + eta seeded: every edge
        gets b and eta b with eta of odd order)
  s_e letters per edge, s_in letters per input per move (inputs offered to moves with codepth >= cmin)
  dlist: comma list of layer thresholds to solve, default m..1
Prints per skeleton the residual floors res_d (min over starts) for each d."""
import sys, time, random
from b5eng import *

SHAPES = {
    # 2 in 7 at block 5, survivors of the two-layer killer KT (tshapes.py), up to mirror
    'E5': [(0, 3), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
    'BF': [(0, 2), (1, 2), (2, 3), (3, 4), (4, 5), (4, 6)],
    'E7': [(0, 2), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
    'E12': [(0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6)],
    # depth-tight control (2 in 6, block 5): must fail at d = 4 by Theorem T
    'C6': [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)],
}


def chain(T):
    return [(i, i + 1) for i in range(T - 1)]


def poset_info(T, covers):
    pred = [set() for _ in range(T)]
    for (i, j) in covers:
        pred[j].add(i)
    below = [set() for _ in range(T)]
    for j in range(T):
        for i in pred[j]:
            below[j] |= below[i] | {i}
    L = [[0] * T for _ in range(T)]
    for j in range(T):
        L[j][j] = 1
        for i in range(j):
            L[i][j] = max([L[i][p] + 1 for p in pred[j] if L[i][p]] + [0])
    depth = [max(L[i][j] for i in range(T)) for j in range(T)]
    codepth = [max(L[i][j] for j in range(T)) for i in range(T)]
    return pred, below, depth, codepth


def make_menu(E, kind, rng):
    Lt = Letters(E); m, n = E.m, E.n
    if kind == 'rand':
        return lambda: E.affine() @ E.phase() @ E.Xs[rng.randrange(n)]
    lv = levi_menu(m)
    eta = lv['eta%d' % ((1 << m) - 1)]
    if kind == 'singer':
        mats = [Lt.addr(gl_pow(eta, k, m)) for k in range((1 << m) - 1)]
        return lambda: mats[rng.randrange(len(mats))] @ E.Xs[rng.randrange(n)]
    if kind in ('levi', 'levi1p'):
        gens = []
        for name, P in lv.items():
            o = gl_order(P, m)
            for k in range(o):
                gens.append(Lt.addr(gl_pow(P, k, m)))
        ph = [Lt.sphase(a) for a in range(n)] + [Lt.cz(j, k) for j in range(m) for k in range(j + 1, m)]
        def f():
            M = gens[rng.randrange(len(gens))]
            if rng.random() < 0.5:
                M = M @ ph[rng.randrange(len(ph))]
            return M @ E.Xs[rng.randrange(n)]
        f.eta = Lt.addr(eta); f.lv = {k: Lt.addr(v) for k, v in lv.items()}
        return f
    raise ValueError(kind)


def leads(E, P, gi):
    """per array l: the term of generator gi with array l closest to the kernel frame (min dK)"""
    out = {}
    for (l, g) in P.gens[gi]:
        d = E.dK(g)
        if l not in out or d < out[l][0]:
            out[l] = (d, g)
    return out


def step_ok(E, b, g):
    h = E.intern(E.snap(E.D @ b @ E.MATS[g]))
    return E.dK(h) == E.dK(g) - 1


def pick_letters(E, menu, gl, rng, want, tries=300):
    """gl: list of lead group elements (one per array).  Prefer letters progressing ALL of them (joint geodesic letter,
    Lemma L), then letters progressing at least one, cycling over the arrays."""
    joint, single = [], {i: [] for i in range(len(gl))}
    for _ in range(tries):
        b = menu()
        ok = [step_ok(E, b, g) for g in gl]
        if all(ok):
            joint.append(b)
            if len(joint) >= want:
                break
        else:
            for i, o in enumerate(ok):
                if o and len(single[i]) < want:
                    single[i].append(b)
    out = joint[:want]
    i = 0
    while len(out) < want and any(single.values()):
        if single[i % len(gl)]:
            out.append(single[i % len(gl)].pop())
        i += 1
    return out, len(joint)


def build(E, W, T, covers, menu, s_e, s_in, cmin, rng, kind, prog=True):
    pred, below, depth, codepth = poset_info(T, covers)
    P = Prog(E, W)
    gidx = []
    njoint = 0
    for k in range(T):
        summ = []
        for i in sorted(pred[k]):
            L = leads(E, P, gidx[i])
            gl = [L[l][1] for l in sorted(L)]
            bs, nj = pick_letters(E, menu, gl, rng, s_e) if prog else ([], 0)
            njoint += min(nj, 1)
            while len(bs) < s_e:
                bs.append(menu())
            for s, b in enumerate(bs):
                summ.append((gidx[i], b))
                if kind == 'levi1p' and s == 0:
                    etas = list(menu.lv.values())
                    summ.append((gidx[i], etas[rng.randrange(len(etas))] @ b))
        if codepth[k] >= cmin or not pred[k]:
            for l in range(W):
                n_in = s_in if pred[k] else max(1, s_in)
                bs, _ = pick_letters(E, menu, [E.ID], rng, n_in) if prog else ([], 0)
                while len(bs) < n_in:
                    bs.append(menu())
                for b in bs:
                    summ.append((l, b))
        gidx.append(P.add_move(summ, desc=(k, sorted(pred[k]))))
    P.njoint = njoint
    return P, depth, codepth


def main():
    a = sys.argv
    shape, W, m, kind = a[1], int(a[2]), int(a[3]), a[4]
    s_e, s_in, starts, seed, nskel = int(a[5]), int(a[6]), int(a[7]), int(a[8]), int(a[9])
    cap = int(a[10]) if len(a) > 10 else 2000
    dl = [int(x) for x in a[11].split(',')] if len(a) > 11 else list(range(m, 0, -1))
    cmin = int(a[12]) if len(a) > 12 else 1
    if shape.startswith('chain'):
        T = int(shape[5:]); covers = chain(T)
    else:
        covers = SHAPES[shape]; T = 1 + max(j for (_, j) in covers)
    rng = random.Random(seed)
    t0 = time.time()
    allres = []
    for s in range(nskel):
        E = Eng(m, seed=seed * 1000 + s)          # fresh engine per skeleton (bounded memory)
        menu = make_menu(E, kind, rng)
        P, depth, codepth = build(E, W, T, covers, menu, s_e, s_in, cmin, rng, kind)
        nterms = [len(g) for g in P.gens]
        if not P.setup_final(cap, dmin=min(dl)):
            print("skel %d: refused (%s) (terms %s)" % (s, getattr(P, "why", "?"), nterms[W:]), flush=True); continue
        fi = P.free_idx(); nf = len(fi)
        dims = [(f[0], f[1]) for f in P.fin]
        res = {}
        xbest = {}
        for d in dl:
            best = np.inf
            for st in range(starts):
                x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
                if nf == 0:
                    r, _ = P.vp(np.ones(P.nal, dtype=complex), d, [], want_J=False); c = float(np.sum(np.abs(r) ** 2))
                    x = x0
                else:
                    x, c = P.lm(d, x0, fi)
                if c < best:
                    best = c; xbest[d] = x
                if best < 1e-20:
                    break
            res[d] = best
            if best > 1e-8:
                break          # deeper layers cannot pass if this one does not
        import resource
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1 << 20 if sys.platform == 'darwin' else 1 << 10)
        print("skel %d: shape %s menu %s s_e %d s_in %d unknowns %d joint %d terms %s closure %s | %s  (%.0fs, %d MB, %d elts)" % (
            s, shape, kind, s_e, s_in, nf, P.njoint, nterms[W:], dims,
            " ".join("d%d=%.2e" % (d, res[d]) for d in sorted(res, reverse=True)), time.time() - t0, rss, len(E.MATS)), flush=True)
        allres.append(res)
        dmin = min([d for d in res if res[d] < 1e-12], default=None)
        if dmin is not None and dmin <= 3:
            np.save("sol_%s_%s_%d_%d_d%d.npy" % (shape, kind, seed, s, dmin), xbest[dmin])
    # summary: deepest layer passed per skeleton
    deep = [min([d for d in r if r[d] < 1e-12], default=m + 1) for r in allres]
    from collections import Counter
    print("SUMMARY shape %s W %d m %d menu %s: deepest passed layer counts %s" % (shape, W, m, kind, sorted(Counter(deep).items())), flush=True)


if __name__ == "__main__":
    main()
