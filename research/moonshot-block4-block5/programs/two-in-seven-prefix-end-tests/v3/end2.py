"""end2.py -- o3-build: EXACT solver for a shared two-move end with group-algebra coefficients.
Given prefix generators (inputs first), paid elements E6, E7 and coefficient patterns alpha[o][j] (lists of
(rational, free letter matrix)), decide whether contents c6, c7 (Q[B']-combinations of prefix generators, SHARED by
both outputs) and per-output combinations sigma_o exist with
     K x_o = sigma_o + sum_j sum_{(a, lam) in alpha[o][j]} a * lam * E_j * c_j        (o = 0, 1).
This is linear in (c6, c7, sigma_0, sigma_1) once the patterns are fixed.  Patterns used:
  PRIV  [[1],[0]] / [[0],[1]]                      two private ends (= free_last_test for each array)
  HALF  [[1, 1], [1, -1]]                          +-1 halving end (E6 = E7: change of basis of PRIV)
  O3    [[(1+s-s^2), (1-s+s^2)], [(1-s+s^2), -(1-s+s^2)]]   order-3 halving end: u6 ~ K(x1+x2), u7 ~ K(x1 - s x2);
        NOT a change of basis of PRIV when s does not normalise E (controls.py C3).
Rows = (output, array, group element); free terms are dropped (cancelled by free images of the inputs in sigma_o).
Solved mod p (bm.rref, C) and then exactly over Q.  Returns ("yes", solution) / ("no", None) / ("cap", None)."""
import numpy as np
from fractions import Fraction
import bm as BM

def joint_closure(G, gens, Es, alpha, cap, outs=(0, 1)):
    index = {}
    for i, g in enumerate(gens):
        for (ll, gid) in g:
            index.setdefault((G.frm(gid), ll), []).append((i, gid))
    # pre-multipliers per (o, j): list of (coef, M = lam E_j, Minv)
    pres = {}
    for o in (0, 1):
        for j, E in enumerate(Es):
            lst = []
            for (a, lam) in alpha[o][j]:
                M = G.mul(lam, E)
                lst.append((Fraction(a), M, np.conj(M).T))
            pres[(o, j)] = lst
    rows = {}; cols = {}; ent = []; queue = []
    def addrow(t):
        r = rows.get(t)
        if r is None:
            r = len(rows); rows[t] = r; queue.append(t)
        return r
    rowcap = 3 * cap
    for o in outs:
        addrow((o, o, G.KID))
    while queue:
        if len(rows) > rowcap:
            return None
        (o, lt, gt) = queue.pop()
        Mg = G.MATS[gt]
        # sigma_o columns
        for (i, gterm) in index.get((G.frm(gt), lt), ()):
            bmat = G.mul(Mg, np.conj(G.MATS[gterm]).T)
            if not G.is_monomial(bmat):
                continue
            b = G.intern(bmat)
            ck = ('S', o, i, b)
            if ck in cols:
                continue
            c = len(cols); cols[ck] = c
            if c >= cap:
                return None
            for (l2, g2), co in gens[i].items():
                g3 = G.intern(G.mul(bmat, G.MATS[g2]))
                if G.is_free(g3):
                    continue
                ent.append((addrow((o, l2, g3)), c, co))
        # shared content columns
        for j in range(len(Es)):
            for (a, M, Mi) in pres[(o, j)]:
                hm = G.mul(Mi, Mg)
                fk = G.frame_key(hm)
                for (i, gterm) in index.get((fk, lt), ()):
                    bmat = G.mul(hm, np.conj(G.MATS[gterm]).T)
                    if not G.is_monomial(bmat):
                        continue
                    b = G.intern(bmat)
                    ck = ('E', j, i, b)
                    if ck in cols:
                        continue
                    c = len(cols); cols[ck] = c
                    if c >= cap:
                        return None
                    for oo in outs:
                        for (a2, M2, _) in pres[(oo, j)]:
                            pre = G.mul(M2, bmat)
                            for (l2, g2), co in gens[i].items():
                                g3 = G.intern(G.mul(pre, G.MATS[g2]))
                                if G.is_free(g3):
                                    continue
                                ent.append((addrow((oo, l2, g3)), c, co * a2))
    return rows, cols, ent

def solve(nr, nc, ent, exact=True, ntgt=2):
    P = BM.P
    if nc == 0:
        return None
    A = np.zeros((nr, nc + 1))
    for (r, c, co) in ent:
        co = Fraction(co)
        A[r, c] = (A[r, c] + (co.numerator % P) * pow(co.denominator % P, P - 2, P)) % P
    for t in range(ntgt):
        A[t, nc] = 1
    Rr, piv = BM.rref(A, full=False)
    if nc in piv:
        return None
    if not exact:
        return True
    # exact over Q restricted to rows/cols that matter (dense Fractions; sizes here are a few hundred)
    M = [[Fraction(0)] * (nc + 1) for _ in range(nr)]
    for (r, c, co) in ent:
        M[r][c] += Fraction(co)
    for t in range(ntgt):
        M[t][nc] = Fraction(1)
    pr = 0; pcols = []
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
        pcols.append(c); pr += 1
        if pr == nr:
            break
    if any(M[r][nc] for r in range(pr, nr)):
        return None
    x = [Fraction(0)] * nc
    for r, c in enumerate(pcols):
        x[c] = M[r][nc]
    return x

def end2(G, gens, Es, alpha, cap=3000, exact=True, outs=(0, 1)):
    cl = joint_closure(G, gens, Es, alpha, cap, outs)
    if cl is None:
        return "cap", None, 0
    rows, cols, ent = cl
    x = solve(len(rows), len(cols), ent, exact, ntgt=len(outs))
    if x is None:
        return "no", None, len(cols)
    if x is True:
        return "yes", None, len(cols)
    inv = {c: k for k, c in cols.items()}
    return "yes", [(inv[c], x[c]) for c in range(len(x)) if x[c]], len(cols)

def check_solution(G, gens, Es, alpha, sol):
    """formal re-check: rebuild sum per output and compare with K x_o up to free terms"""
    acc = [{}, {}]
    for (key, x) in sol:
        if key[0] == 'S':
            _, o, i, b = key
            for (l2, g2), co in gens[i].items():
                t = (l2, G.intern(G.mul(G.MATS[b], G.MATS[g2])))
                acc[o][t] = acc[o].get(t, 0) + x * co
        else:
            _, j, i, b = key
            for o in (0, 1):
                for (a, lam) in alpha[o][j]:
                    pre = G.mul(lam, Es[j], G.MATS[b])
                    for (l2, g2), co in gens[i].items():
                        t = (l2, G.intern(G.mul(pre, G.MATS[g2])))
                        acc[o][t] = acc[o].get(t, 0) + x * co * Fraction(a)
    ok = True
    for o in (0, 1):
        d = {t: v for t, v in acc[o].items() if v != 0}
        if d.pop((o, G.KID), None) != 1:
            ok = False
        if not all(G.is_free(g) for (_, g) in d):
            ok = False
    return ok

def patterns(G, s):
    I = G.I
    s2 = G.mul(s, s)
    one = [(1, I)]
    A = [(1, I), (1, s), (-1, s2)]        # 1 + s - s^2
    B = [(1, I), (-1, s), (1, s2)]        # 1 - s + s^2
    Bm = [(-1, I), (1, s), (-1, s2)]
    return {
        "PRIV": [[one, []], [[], one]],
        "HALF": [[one, one], [one, [(-1, I)]]],
        "O3": [[A, B], [B, Bm]],
    }

def raw_check_end2(w, prog, gens, Es, alpha, sol, fs=(1, 2, 3), seed=7):
    """replay the prefix on random complex arrays at each f, form u6 = E6 c6, u7 = E7 c7 and the outputs
    sum_j alpha[o][j] (lam E_j c_j) + sigma_o, subtract the formal free leftover (free images of the inputs),
    compare with K^(x)f x_o.  Returns {f: max relative error}."""
    from o3lib import apply_tensor, N_
    G = w.G
    rng = np.random.default_rng(seed)
    out = {}
    # formal leftover per output
    lefts = [{}, {}]
    for (key, x) in sol:
        if key[0] == 'S':
            _, o, i, b = key
            for (l2, g2), co in gens[i].items():
                t = (l2, G.intern(G.mul(G.MATS[b], G.MATS[g2])))
                lefts[o][t] = lefts[o].get(t, 0) + x * co
        else:
            _, j, i, b = key
            for o in (0, 1):
                for (a, lam) in alpha[o][j]:
                    pre = G.mul(lam, Es[j], G.MATS[b])
                    for (l2, g2), co in gens[i].items():
                        t = (l2, G.intern(G.mul(pre, G.MATS[g2])))
                        lefts[o][t] = lefts[o].get(t, 0) + x * co * Fraction(a)
    for o in (0, 1):
        lefts[o][(o, G.KID)] = lefts[o].get((o, G.KID), 0) - 1
        lefts[o] = {t: v for t, v in lefts[o].items() if v != 0}
        assert all(G.is_free(g) for (_, g) in lefts[o]), "non-free leftover"
    for f in fs:
        Nf = N_ ** f
        xs = [rng.normal(size=Nf) + 1j * rng.normal(size=Nf) for _ in range(2)]
        vecs = list(xs)
        for mv in prog:
            c = np.zeros(Nf, dtype=complex)
            for (nu, de, a, j) in mv["cs"]:
                v = vecs[j] if a is None else apply_tensor(G.mono(*a), vecs[j], f)
                c = c + (nu / de) * v
            vecs.append(apply_tensor(G.mul(w.D1, G.mono(*mv["b"])), c, f))
        cs = [np.zeros(Nf, dtype=complex) for _ in Es]
        sig = [np.zeros(Nf, dtype=complex), np.zeros(Nf, dtype=complex)]
        for (key, x) in sol:
            if key[0] == 'S':
                _, o, i, b = key
                sig[o] += float(x) * apply_tensor(G.MATS[b], vecs[i], f)
            else:
                _, j, i, b = key
                cs[j] += float(x) * apply_tensor(G.MATS[b], vecs[i], f)
        us = [apply_tensor(E, c, f) for E, c in zip(Es, cs)]
        err = 0.0
        for o in (0, 1):
            y = sig[o].copy()
            for j in range(len(Es)):
                for (a, lam) in alpha[o][j]:
                    y += float(a) * apply_tensor(lam, us[j], f)
            for (l2, g), v in lefts[o].items():
                y -= float(v) * apply_tensor(G.MATS[g], xs[l2], f)
            tgt = apply_tensor(G.K, xs[o], f)
            err = max(err, float(np.abs(y - tgt).max() / np.abs(tgt).max()))
        out[f] = err
    return out
