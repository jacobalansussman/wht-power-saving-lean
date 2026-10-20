"""joint2.py -- ms-fano: EXACT end-gadget completion with TWO SHARED extra moves.
Given a prefix (generators), ask whether there are contents c1, c2 (Q[B']-combinations of the prefix generators and
raw inputs, solved) such that, with the two extra generators E1 c1 and E2 c2,
    K x1 = span + E1 c1 + E2 c2      and      K x2 = span' + a1 E1 c1 + a2 E2 c2
(span, span' = Q[B']-combinations of the prefix); more generally output o uses C[o][j] E_j c_j with the coefficient
matrix C one of [[1,0],[0,1]] (two private-content last moves, as fano7), [[1,1],[1,-1]] (the halving end gadget:
E1 c1 = (K x1 + K x2)/2-type), [[1,1],[0,1]] (one move shared by both outputs).  For fixed E1, E2, a1, a2 this is ONE linear system (columns of c1, c2 shared by both outputs); it is built by
a joint closure over frames (as tb4-beam bm.Game.closure, two outputs, two paid columns) and solved mod p, then exactly
over Q.  Moves of the word = prefix + 2.  Block 3: a 3-move prefix can never pass (2 in 5); block 4: 5-move prefixes
give 2 in 7.
usage: python3 joint2.py m kind k n seed [nE] [out.json]   (kind as in fano7.py: PL FS FN FE)"""
import sys, time, json, random
from fractions import Fraction
from fano7 import *
from stag2 import junk_E
import bm as BM

def joint_closure(G, gens, Es, alphas, cap):
    index = {}
    for i, g in enumerate(gens):
        for (ll, gid) in g:
            index.setdefault((G.frm(gid), ll), []).append((i, gid))
    Einvs = [np.conj(E).T for E in Es]
    def freeish(g3):
        if G.is_free(g3):
            return True
        return any(G.is_monomial(G.mul(Ei, G.MATS[g3])) for Ei in Einvs)
    t0 = (0, 0, G.KID); t1 = (1, 1, G.KID)
    rows = {t0: 0, t1: 1}; cols = {}; ent = []; queue = [t0, t1]
    def addrow(t):
        r = rows.get(t)
        if r is None:
            r = len(rows); rows[t] = r; queue.append(t)
        return r
    while queue:
        (o, lt, gt) = queue.pop()
        cands = []
        for (i, gterm) in index.get((G.frm(gt), lt), ()):
            cands.append(((o, 'S', i), G.mul(G.MATS[gt], np.conj(G.MATS[gterm]).T), None))
        for j, Ei in enumerate(Einvs):
            if alphas[o][j] == 0:
                continue
            fk = G.frame_key(G.mul(Ei, G.MATS[gt]))
            for (i, gterm) in index.get((fk, lt), ()):
                cands.append(((j, 'E', i), G.mul(Ei, G.MATS[gt], np.conj(G.MATS[gterm]).T), j))
        for (key, bmat, j) in cands:
            if not G.is_monomial(bmat):
                continue
            b = G.intern(bmat)
            ck = key + (b,)
            if ck in cols:
                continue
            c = len(cols); cols[ck] = c
            if c >= cap:
                return None
            if j is None:
                for (l2, g2), co in gens[key[2]].items():
                    g3 = G.intern(G.mul(bmat, G.MATS[g2]))
                    if freeish(g3):
                        continue
                    ent.append((addrow((key[0], l2, g3)), c, co))
            else:
                pre = G.mul(Es[j], bmat)
                for (l2, g2), co in gens[key[2]].items():
                    g3 = G.intern(G.mul(pre, G.MATS[g2]))
                    if freeish(g3):
                        continue
                    for oo in (0, 1):
                        if alphas[oo][j] != 0:
                            ent.append((addrow((oo, l2, g3)), c, co * alphas[oo][j]))
    return rows, cols, ent

def solve(nr, nc, ent, exact=True):
    P = BM.P
    A = np.zeros((nr, nc + 1))
    for (r, c, co) in ent:
        co = Fraction(co)
        A[r, c] = (A[r, c] + (co.numerator % P) * pow(co.denominator % P, P - 2, P)) % P
    A[0, nc] = 1; A[1, nc] = 1
    Rr, piv = BM.rref(A, full=False)
    if nc in piv:
        return None
    if not exact:
        return True
    # exact over Q (sparse-ish Gauss on Fractions)
    M = [[Fraction(0)] * (nc + 1) for _ in range(nr)]
    for (r, c, co) in ent:
        M[r][c] += Fraction(co)
    M[0][nc] = Fraction(1); M[1][nc] = Fraction(1)
    pr = 0
    for c in range(nc):
        p_ = next((r for r in range(pr, nr) if M[r][c]), None)
        if p_ is None:
            continue
        M[pr], M[p_] = M[p_], M[pr]
        inv = 1 / M[pr][c]
        M[pr] = [x * inv for x in M[pr]]
        for r in range(nr):
            if r != pr and M[r][c]:
                f = M[r][c]; M[r] = [x - f * y for x, y in zip(M[r], M[pr])]
        pr += 1
    if any(M[r][nc] for r in range(pr, nr)):
        return None
    return True

def E_pool(ft, gens, nE):
    """paid elements E = K g^-1 for every prefix term g one step from the kernel frame (so E g = K exactly), for
    either array, plus junk-cancelling ones (stag2.junk_E); at most 2 nE, preferring distinct frames"""
    G = ft.G
    out = {}
    for gen in gens:
        for (l, gid) in gen:
            if G.dist(gid) != 1:
                continue
            E = G.mul(G.K, np.conj(G.MATS[gid]).T)
            if not (G.rowsupp(E) == 2).all():
                continue
            out.setdefault(G.key(E), E)
    for l in range(2):
        for E in junk_E(ft, gens, l, nE):
            out.setdefault(G.key(E), E)
    vals = list(out.values())
    random.Random(len(vals)).shuffle(vals)
    return vals[:2 * nE]

def _deadline():
    import os
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'DEADLINE')
    try:
        return float(open(f).read().split()[0])
    except Exception:
        return None

def _budget():
    import os
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BUDGET')
    try:
        return float(open(f).read().split()[0])
    except Exception:
        return None

def main():
    m, kind, k, n, seed = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    nE = int(sys.argv[6]) if len(sys.argv) > 6 else 4
    out = sys.argv[7] if len(sys.argv) > 7 else None
    t0 = time.time(); rng = random.Random(seed)
    ft = setup(m, seed)
    st = {"prefixes": 0, "fail_build": 0, "systems": 0, "cap": 0, "pass": 0, "pass_prefixes": 0}
    hits = []
    dl = _deadline()
    for it in range(n):
        bu = _budget()
        if (dl is not None and time.time() > dl) or (bu is not None and time.time() - t0 > bu):
            st['stopped_at_deadline_after'] = it; break
        if it and it % 10 == 0:
            print('progress', it, st, flush=True)
        if len(ft.G.MATS) > 40000:
            ft = setup(m, seed + 1000 + it)
        gens = prefix(ft, kind, 2, k, rng)
        if gens is None:
            st["fail_build"] += 1; continue
        st["prefixes"] += 1
        Es = E_pool(ft, gens, nE)
        found = False
        for a in range(len(Es)):
            for b in range(len(Es)):
                for alphas in (((1, 0), (0, 1)), ((1, 1), (1, -1)), ((1, 1), (0, 1))):
                    st["systems"] += 1
                    cl = joint_closure(ft.G, gens, [Es[a], Es[b]], alphas, 1200)
                    if cl is None:
                        st["cap"] += 1; continue
                    rows, cols, ent = cl
                    if solve(len(rows), len(cols), ent):
                        st["pass"] += 1; found = True; hits.append((it, a, b, alphas))
                        break
                if found: break
            if found: break
        st["pass_prefixes"] += found
    st["seconds"] = round(time.time() - t0, 1)
    print(f"joint2 m={m} {kind} k={k} n={n} seed={seed} nE={nE}: {st}", flush=True)
    if out:
        json.dump({"args": sys.argv[1:], "stats": st, "hits": hits}, open(out, "w"))

if __name__ == "__main__":
    main()
