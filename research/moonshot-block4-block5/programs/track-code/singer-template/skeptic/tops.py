"""tops.py (mw2-singer4): the SEPARATION-MATRIX end test, all T at once.

Given a prefix (generators of M - W moves) and a paid element E (or None), compute the space
    V_E = { t in C[H]^W : sum_o t_o K x_o  lies in  E * C[B']{prefix, inputs} + C[B']{prefix} }
for a finite group H of free letters (default H = <eta>, the Singer cycle of F_2^4/<x> lifted to the Levi; H = {1} is the
scalar case).  V_E is a linear subspace of C^{|H| W}: one closure over frames from the |H| W target rows (o, hK), one
echelon reduction mod p of [A | -P] (rows with a pivot in a t-column are the constraints on t).
A unit "W in M" of Singer-separation form = W tops E_1..E_W (repetition allowed) and t^(j) in V_{E_j} whose W x W matrix
T over C[H] is invertible (K x = T^-1 (tops + span), T^-1 in Q[H] subset Q[B']).  For H abelian (Z7) T is invertible iff
det T(chi) != 0 for every character chi; a generic choice t^(j) in V_{E_j} achieves this iff for every chi the spaces
chi(V_{E_j}) admit an independent transversal (Rado: dim sum_{j in J} chi(V_{E_j}) >= |J| for all J).  Computed mod
p = 1048573 (p = 1 mod 28: 7th roots of unity and i exist mod p); a mod-p pass is then to be confirmed exactly."""
import itertools, time
import numpy as np
from fractions import Fraction as Fr
import bm as BM
from s4lib import S4, E_pool, fast_frame_key

P = BM.P
OMEGA = next(pow(g, (P - 1) // 7, P) for g in range(2, 1000) if pow(g, (P - 1) // 7, P) != 1)
assert pow(OMEGA, 7, P) == 1

def frac_p(co):
    co = Fr(co)
    return (co.numerator % P) * pow(co.denominator % P, P - 2, P) % P

def target_space(S, gens, W, E=None, Hids=None, cap=2000, rowcap=6000):
    """returns (basis of V_E as array (dim, |H| W) mod p, stats) or (None, stats) on cap"""
    G = S.G
    if Hids is None:
        Hids = S.etap_id
    nh = len(Hids)
    index = {}
    for i, g in enumerate(gens):
        for (l, gid) in g:
            index.setdefault((S.fk(gid), l), []).append((i, gid))
    Eid = None if E is None else G.intern(E)
    Einv = None if E is None else S.inv(Eid)
    rows = {}; cols = {}; ent = []; queue = []
    def addrow(t):
        r = rows.get(t)
        if r is None:
            if len(rows) >= rowcap:
                raise OverflowError
            r = len(rows); rows[t] = r; queue.append(t)
        return r
    targets = []
    try:
        for o in range(W):
            for h in Hids:
                targets.append(addrow((o, S.mulid(h, G.KID))))
        while queue:
            (l, g) = queue.pop()
            if len(G.MATS) > 150000:
                raise OverflowError
            cands = []
            for (i, gt) in index.get((S.fk(g), l), ()):
                cands.append((('S', i), S.mono_id(g, S.inv(gt))))
            if Eid is not None:
                hM = S.prod(Einv, g)
                for (i, gt) in index.get((fast_frame_key(hM), l), ()):
                    cands.append((('E', i), S.mono_id(Einv, g, S.inv(gt))))
            for (key, b) in cands:
                if b is None:
                    continue
                ck = key + (b,)
                if ck in cols:
                    continue
                c = len(cols); cols[ck] = c
                if c >= cap:
                    return None, {'rows': len(rows), 'cols': len(cols), 'cap': True}
                pre = b if key[0] == 'S' else S.mulid(Eid, b)
                for (l2, g2), co in gens[key[1]].items():
                    g3 = S.mulid(pre, g2)
                    if S.is_mono(g3):
                        continue
                    ent.append((addrow((l2, g3)), c, co))
    except OverflowError:
        return None, {'rows': len(rows), 'cols': len(cols), 'cap': True}
    nr, nc = len(rows), len(cols); q = nh * W
    A = np.zeros((nr, nc + q))
    for (r, c, co) in ent:
        A[r, c] = (A[r, c] + frac_p(co)) % P
    for idx, r in enumerate(targets):
        A[r, nc + idx] = (A[r, nc + idx] - 1) % P
    Rr, piv = BM.rref(A, full=False)
    cons = [i for i, pc in enumerate(piv) if pc >= nc]
    C = Rr[cons][:, nc:] if cons else np.zeros((0, q))
    if C.shape[0] == 0:
        V = np.eye(q)
    else:
        V = BM.nullspace(C).T          # rows = basis vectors of V
    return V % P, {'rows': nr, 'cols': nc, 'dimV': V.shape[0], 'q': q}

def char_image(V, W, nh, s):
    """chi_s(V): t (o, k) -> sum_k t_{o,k} omega^{s k}; V rows (dim, nh W) with index o * nh + k"""
    if V is None or V.shape[0] == 0:
        return np.zeros((0, W))
    if nh == 1:
        return V.copy()
    w = np.array([pow(OMEGA, (s * k) % 7, P) for k in range(nh)], dtype=np.float64)
    out = np.zeros((V.shape[0], W))
    for o in range(W):
        out[:, o] = BM.mmp(V[:, o * nh:(o + 1) * nh], w[:, None])[:, 0]
    return out

def rank_p(M):
    if M.shape[0] == 0:
        return 0
    return BM.rank(M % P)

def rado_ok(images, W):
    """images: list of W arrays (rows span chi(V_j)); independent transversal exists?"""
    for r in range(1, W + 1):
        for J in itertools.combinations(range(W), r):
            M = np.concatenate([images[j] for j in J], axis=0)
            if rank_p(M) < r:
                return False
    return True

def best_tops(spaces, W, nh):
    """spaces: dict name -> V.  Find a multiset of W names passing Rado for every character. Returns (names, per-char
    diagnostics) or (None, diagnostics): diagnostic = for each character, the max rank of the sum of all images."""
    names = [k for k, V in spaces.items() if V is not None and V.shape[0] > 0]
    nchar = 7 if nh == 7 else 1
    imgs = {k: [char_image(spaces[k], W, nh, s) for s in range(nchar)] for k in names}
    diag = []
    for s in range(nchar):
        tot = np.concatenate([imgs[k][s] for k in names], axis=0) if names else np.zeros((0, W))
        diag.append(rank_p(tot))
    if min(diag) < W:
        return None, diag
    for combo in itertools.combinations_with_replacement(names, W):
        if all(rado_ok([imgs[k][s] for k in combo], W) for s in range(nchar)):
            return combo, diag
    return None, diag

def target_space_full(S, gens, W, E, cap=2000, rowcap=6000):
    """ALL separation matrices at once (T over C[B'], not only C[<eta>]): the space of K-frame parts t (coefficients on
    the K-frame elements b K x_o reached) such that E c + span is supported on the K-frame (c, span solved).  Seeds:
    top columns b * gen_i with E b g on the K-frame.  Returns (augmentation image of V as array (dim, W) mod p, stats):
    the D2 test (top augmentations) needs W independent rows over all tops; this is the exact-closure relaxation gate."""
    G = S.G
    index = {}
    for i, g in enumerate(gens):
        for (l, gid) in g:
            index.setdefault((S.fk(gid), l), []).append((i, gid))
    Eid = G.intern(E); Einv = S.inv(Eid)
    Kfk = S.fk(G.KID)
    rows = {}; cols = {}; ent = []; queue = []
    def addrow(t):
        r = rows.get(t)
        if r is None:
            if len(rows) >= rowcap:
                raise OverflowError
            r = len(rows); rows[t] = r
            if S.fk(t[1]) != Kfk:
                queue.append(t)
        return r
    def addcol(key, b):
        ck = key + (b,)
        if ck in cols:
            return
        c = len(cols); cols[ck] = c
        if c >= cap:
            raise OverflowError
        pre = b if key[0] == 'S' else S.mulid(Eid, b)
        for (l2, g2), co in gens[key[1]].items():
            g3 = S.mulid(pre, g2)
            if S.is_mono(g3):
                continue
            ent.append((addrow((l2, g3)), c, co))
    try:
        # seeds: E b g = K exactly for prefix terms g with K-frame-completing E
        Kinv_target = S.mulid(Einv, G.KID)          # E^-1 K
        fkt = S.fk(Kinv_target)
        for l in range(W):
            for (i, gt) in index.get((fkt, l), ()):
                b = S.mono_id(Kinv_target, S.inv(gt))
                if b is not None:
                    addcol(('E', i), b)
        while queue:
            if len(G.MATS) > 150000:
                raise OverflowError
            (l, g) = queue.pop()
            for (i, gt) in index.get((S.fk(g), l), ()):
                b = S.mono_id(g, S.inv(gt))
                if b is not None:
                    addcol(('S', i), b)
            hM = S.prod(Einv, g)
            for (i, gt) in index.get((fast_frame_key(hM), l), ()):
                b = S.mono_id(Einv, g, S.inv(gt))
                if b is not None:
                    addcol(('E', i), b)
    except OverflowError:
        return None, {'rows': len(rows), 'cols': len(cols), 'cap': True}
    nr, nc = len(rows), len(cols)
    if nc == 0:
        return np.zeros((0, W)), {'rows': nr, 'cols': 0}
    Krows = [(t, r) for t, r in rows.items() if S.fk(t[1]) == Kfk]
    # unknowns x (columns); K-frame rows are free outputs t = A_K x; other rows: A_N x = 0.
    nonK = [r for t, r in rows.items() if S.fk(t[1]) != Kfk]
    A = np.zeros((nr, nc))
    for (r, c, co) in ent:
        A[r, c] = (A[r, c] + frac_p(co)) % P
    AN = A[nonK] if nonK else np.zeros((0, nc))
    X = BM.nullspace(AN) if AN.shape[0] else np.eye(nc)      # columns: kernel basis
    if X.shape[1] == 0:
        return np.zeros((0, W)), {'rows': nr, 'cols': nc, 'ker': 0}
    # augmentation per array of the K-frame part: aug_o(x) = sum over K rows of array o of (A x)
    Aug = np.zeros((W, nc))
    for (t, r) in Krows:
        Aug[t[0]] = (Aug[t[0]] + A[r]) % P
    img = BM.mmp(Aug, X).T % P          # (ker dim, W)
    return img, {'rows': nr, 'cols': nc, 'ker': int(X.shape[1]), 'Krows': len(Krows), 'augrank': rank_p(img)}
