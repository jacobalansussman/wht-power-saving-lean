"""lsearch2.py (mw2-b5beyond): DESIGNED skeletons + layered numerical solve (block m, W arrays).

Skeleton (any dependency shape): along every cover edge the main letter is a geodesic step for the leads of the arrays
present (joint letter if one exists, else a single-array step); top moves (depth >= m) are ASSIGNED to arrays
round-robin and get a SPLIT letter: progress for their array's lead, non-progress for the other arrays' leads (so the
top system is diagonal, the Theorem-T square).  Then s_side extra summands per move (random earlier generator or input,
letter from a structured menu: odd-order Levi x Clifford phase x Pauli, or the 1 + eta separation of the main letter).
Solve: variable-projection LM on the residual restricted to layers >= d, d = m, m-1, ... (stop at first failure).
usage: lsearch2.py SHAPE W m side_menu s_side starts seed nskel [sep]   (sep=1: add eta*b next to every main letter)"""
import sys, time, random, resource
from collections import Counter
from b5eng import *
from lsearch import SHAPES, chain, poset_info, make_menu


def leads(E, P, gi):
    out = {}
    for (l, g) in P.gens[gi]:
        d = E.dK(g)
        if l not in out or d < out[l][0]:
            out[l] = (d, g)
    return out


def dstep(E, b, g):
    h = E.intern(E.snap(E.D @ b @ E.MATS[g]))
    return E.dK(g) - E.dK(h)          # +1 progress, 0 level, -1 back


def rand_free(E, rng):
    return E.affine() @ E.phase() @ E.Xs[rng.randrange(E.n)]


def main_letter(E, rng, gl, assign=None, tries=600):
    """gl: dict l -> lead gid.  assign None: prefer joint progress for all; else progress for gl[assign] and
    NO progress for the others."""
    best = None
    for _ in range(tries):
        b = rand_free(E, rng)
        st = {l: dstep(E, b, g) for l, g in gl.items()}
        if assign is None:
            sc = sum(1 for v in st.values() if v == 1)
            if sc == len(gl):
                return b, 'joint'
            if best is None or sc > best[0]:
                best = (sc, b)
        else:
            if st.get(assign) == 1 and all(v < 1 for l, v in st.items() if l != assign):
                return b, 'split'
            if st.get(assign) == 1 and best is None:
                best = (1, b)
    return (best[1], 'partial') if best else (rand_free(E, rng), 'none')


def build(E, W, T, covers, side, s_side, rng, sep, m):
    pred, below, depth, codepth = poset_info(T, covers)
    P = Prog(E, W); gidx = []; tags = []
    tops = [k for k in range(T) if depth[k] >= m]
    starts = [k for k in range(T) if not pred[k]]
    inp_of = {}
    # inputs: one per array at moves of codepth >= m (round robin), others get none in the main skeleton
    S = [k for k in range(T) if codepth[k] >= m]
    for j, k in enumerate(S):
        inp_of.setdefault(k, []).append(j % W)
    for k in range(T):
        summ = []
        for l in inp_of.get(k, []):
            b, tg = main_letter(E, rng, {'x': E.ID}, None)
            summ.append((l, b)); tags.append(tg)
        for i in sorted(pred[k]):
            L = leads(E, P, gidx[i]); gl = {l: L[l][1] for l in L}
            assign = None
            if k in tops and len(gl) > 1:
                assign = sorted(gl)[tops.index(k) % len(gl)]
            b, tg = main_letter(E, rng, gl, assign); tags.append(tg)
            summ.append((gidx[i], b))
            if sep:
                summ.append((gidx[i], side.eta @ b))
        for _ in range(s_side):
            cand = [gidx[i] for i in range(k)] + list(range(W))
            src = cand[rng.randrange(len(cand))]
            summ.append((src, side()))
        if not summ:
            summ.append((rng.randrange(W), rand_free(E, rng)))
        gidx.append(P.add_move(summ, desc=(k, sorted(pred[k]))))
    return P, tags


def main():
    a = sys.argv
    shape, W, m, sidek, s_side, starts, seed, nskel = a[1], int(a[2]), int(a[3]), a[4], int(a[5]), int(a[6]), int(a[7]), int(a[8])
    sep = int(a[9]) if len(a) > 9 else 0
    if shape.startswith('chain'):
        T = int(shape[5:]); covers = chain(T)
    elif shape.startswith('P:'):
        covers = [tuple(map(int, e.split('-'))) for e in shape[2:].split('/')]; T = 1 + max(j for (_, j) in covers)
    else:
        covers = SHAPES[shape]; T = 1 + max(j for (_, j) in covers)
    rng = random.Random(seed)
    t0 = time.time(); deep = []
    for s in range(nskel):
        E = Eng(m, seed=seed * 1000 + s)
        side = make_menu(E, sidek if sidek != 'rand' else 'rand', rng)
        if not hasattr(side, 'eta'):
            side.eta = Letters(E).addr(singer(m))
        P, tags = build(E, W, T, covers, side, s_side, rng, sep, m)
        nterms = [len(g) for g in P.gens][W:]
        if not P.setup_final(cap=800):
            print("skel %d: refused (%s) terms %s" % (s, getattr(P, 'why', '?'), nterms), flush=True); deep.append(m + 2); continue
        fi = P.free_idx(); nf = len(fi); res = {}
        for d in range(m, 0, -1):
            best = np.inf
            for st in range(starts):
                x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
                if nf == 0:
                    r, _ = P.vp(np.ones(P.nal, dtype=complex), d, [], want_J=False); c = float(np.sum(np.abs(r) ** 2))
                else:
                    x, c = P.lm(d, x0, fi, iters=50)
                best = min(best, c)
                if best < 1e-20:
                    break
            res[d] = best
            if best > 1e-10:
                break
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1 << 20 if sys.platform == 'darwin' else 1 << 10)
        dd = min([d for d in res if res[d] < 1e-10], default=m + 1)
        deep.append(dd)
        print("skel %d: %s tags %s unknowns %d terms %s closure %s | %s (%.0fs %d MB)" % (
            s, shape, dict(Counter(tags)), nf, nterms, [(f[0], f[1]) for f in P.fin],
            " ".join("d%d=%.1e" % (d, res[d]) for d in sorted(res, reverse=True)), time.time() - t0, rss), flush=True)
    print("SUMMARY %s W%d m%d side %s s_side %d sep %d: deepest layer passed (m+1 = none, m+2 = refused): %s" % (
        shape, W, m, sidek, s_side, sep, sorted(Counter(deep).items())), flush=True)


if __name__ == "__main__":
    main()
