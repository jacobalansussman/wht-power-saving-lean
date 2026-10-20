"""endsolve.py (w3-b5): EXACT-structure end solver for the 2-in-7 butterfly at block 5 (engine game), with every colliding
injection class offered at once and the end contents left completely free.

Program (mw2-b5beyond PROOF 3b): s_l = D sig_l x_l; t1 = D(a1 s1 + a2 s2 + sum_j alpha_j beta_j x_{l_j}); t2 = D(c2 t1 +
sum_k alpha_k delta_k x_{l_k}); t3 = D(c3 t2); ends e_1, e_2 terminal.  V = C[B]{x1, x2, s1, s2, t1, t2, t3}.
Unit <=> (PROOF section 2, square top, terminal ends) there are c_e in V and kappa in GL_2(C[B]) with D c_e - kappa_e K X in V.
For output l we fix kappa's coefficient on (x_l, K) to 1 and leave every other kappa entry free; the end content c_l and the
V-part are free: unknowns y_{u,b} (D b u in D c), z_{u,b} (b u in V), for u in {t1, t2, t3}; the single-term generators
(x_l, s_l) and kappa make their rows free (exact elimination).  Closure from the target row; rows of layer < d dropped
(layered test).  N(alpha) is LINEAR in the injection coefficients (main letters have coefficient 1), so we minimise
||N(alpha) w - e||^2 by variable projection over alpha only (both outputs share alpha).
Letters are constructed EXACTLY from F_2 symplectic data (sp2.py) and lifted to monomial Cliffords; every lift is checked.
usage: endsolve.py seed lam1 lam2 ninj_per_mech odd_only dlist"""
import sys, time, random, hashlib
import numpy as np
import sp2
from sp2 import (M, NB, MASKX, L0, Fstar, D as DS, K as KS, app, mul as smul, inv as sinv, image, dist, neighbors, rref,
                 span_all, levi, unip, rand_P, order as sorder, IDM)
from mwlib import Eng as Eng0, SC

np.seterr(all="ignore")
ONLY_N2 = False


class Eng(Eng0):
    @staticmethod
    def key(Mx):
        return hashlib.blake2b(np.round(Mx.real * SC).astype(np.int32).tobytes() + np.round(Mx.imag * SC).astype(np.int32).tobytes(), digest_size=16).digest()

    def intern(self, Mx):
        k = self.key(Mx); i = self.IDX.get(k)
        if i is None:
            i = len(self.MATS); self.IDX[k] = i; self.MATS.append(Mx)
        return i


# ------------------------------------------------------------------ symplectic image of an exact matrix
class Symp:
    def __init__(self, E):
        self.E = E; n = E.n
        self.P = []
        for j in range(M):
            self.P.append(E.Xs[1 << j])
        for j in range(M):
            self.P.append(np.diag([(-1.0) ** ((x >> j) & 1) for x in range(n)]).astype(complex))
        self.cache = {}

    def pvec(self, Pm):
        col = Pm[:, 0]; a = int(np.argmax(np.abs(col))); c0 = col[a]
        b = 0
        for j in range(M):
            r = Pm[(1 << j) ^ a, 1 << j] / c0
            if abs(r + 1) < 1e-6:
                b |= 1 << j
            else:
                assert abs(r - 1) < 1e-6, r
        return a | (b << M)

    def of(self, gid):
        v = self.cache.get(gid)
        if v is None:
            G = self.E.MATS[gid]; GH = np.conj(G).T
            v = [self.pvec(G @ P @ GH) for P in self.P]
            self.cache[gid] = v
        return v

    def frame(self, gid):
        """F(g) = gbar^-1 L0"""
        key = ('F', gid)
        v = self.cache.get(key)
        if v is None:
            v = image(sinv(self.of(gid)), L0); self.cache[key] = v
        return v


def lift_P(E, g):
    """monomial Clifford matrix with symplectic image g (g in P = Stab(L0))."""
    n = E.n
    A5 = [g[i] & MASKX for i in range(M)]
    lv = levi(A5)
    u = smul(sinv(lv), g)
    S = [(u[i] >> M) for i in range(M)]
    for i in range(M):
        assert u[i] & MASKX == 1 << i and u[M + i] == 1 << (M + i)
    perm = np.zeros((n, n), complex)
    for x in range(n):
        y = 0
        for j in range(M):
            if (x >> j) & 1:
                y ^= A5[j]
        perm[y, x] = 1
    ph = []
    for x in range(n):
        q = 0
        for j in range(M):
            if (x >> j) & 1:
                q += (S[j] >> j) & 1
                for k in range(j + 1, M):
                    if (x >> k) & 1 and (S[j] >> k) & 1:
                        q += 2
        ph.append(1j ** (q % 4))
    return perm @ np.diag(ph)


def find_P(rng, X, Y, extra=None, tries=400000):
    """random g in P with g X = Y (X, Y frames at the same distance from L0); extra(g) -> bool optional.
    Without `extra`: exact uniform sampling of the coset (fastP.py), no rejection."""
    if extra is None:
        from fastP import fast_find_P
        return fast_find_P(rng, X, Y)
    for t in range(tries):
        g = rand_P(rng)
        if image(g, X) == Y and (extra is None or extra(g)):
            return g
    return None


def opp_graph_T(O):
    """O opposite L0 -> symmetric T with O = {(a, T a)}"""
    S = span_all(O); T = [0] * M
    for v in S:
        a = v & MASKX
        if a and a & (a - 1) == 0:
            T[a.bit_length() - 1] = v >> M
    return T


# ------------------------------------------------------------------ configuration (global picture, NOTE section 1)
def make_config(rng, lam1, lam2):
    M0 = image(DS, L0)
    N1 = rng.choice([N for N in sorted(neighbors(L0)) if dist(M0, N) == 2])
    N2 = rng.choice([N for N in sorted(neighbors(N1)) if dist(L0, N) == 2 and dist(M0, N) == 3])
    anchors = [X for X in sorted(neighbors(M0)) if dist(X, L0) == 2 and dist(X, N1) == 3 and dist(X, N2) == 4]
    nN2 = sorted(neighbors(N2))
    for _ in range(1000):
        P1, Q = rng.sample(anchors, 2)
        A1s = [Y for Y in nN2 if dist(P1, Y) == 5 and dist(Q, Y) == lam1]
        A2s = [Y for Y in nN2 if dist(Q, Y) == 5 and dist(P1, Y) == lam2]
        if A1s and A2s:
            return dict(M0=M0, N1=N1, N2=N2, P1=P1, Q=Q, A1=rng.choice(A1s), A2=rng.choice(A2s), nN2=nN2)
    return None


def build_letters(rng, cf):
    """F_2 letters realising the configuration: c2, c3 (trunk), a1, a2 (merge), s1, s2 (starts, = sigma)."""
    M0, N1, N2 = cf['M0'], cf['N1'], cf['N2']
    c2 = find_P(rng, N1, M0)                                 # c2bar N1 = M0
    Hn = image(DS, image(c2, N2))                            # Dbar c2bar N2, a neighbour of L0
    c3 = find_P(rng, Hn, M0)
    a1 = find_P(rng, M0, image(DS, cf['P1']))                # a1bar M0 = Dbar P1
    a2 = find_P(rng, M0, image(DS, cf['Q']))
    # sigma_l bar F* = Dbar a_l^-1 Dbar A_l  (opposite L0): sigma = unip(T + I)
    out = dict(c2=c2, c3=c3, a1=a1, a2=a2)
    for l, (a, A) in enumerate(((a1, cf['A1']), (a2, cf['A2']))):
        O = image(DS, image(sinv(a), image(DS, A)))
        assert dist(O, L0) == M
        T = opp_graph_T(O)
        S = [T[i] ^ (1 << i) for i in range(M)]
        sg = unip(S)
        assert image(sg, Fstar) == O
        out['s%d' % (l + 1)] = sg
    return out


def pi_bar(L, l):
    """pbar or qbar = Dbar a Dbar sigma"""
    a = L['a%d' % (l + 1)]; s = L['s%d' % (l + 1)]
    return smul(DS, smul(a, smul(DS, s)))


def transports(rng, cf, L, n_per, odd_only):
    """injection letters for every collision mechanism with the two leftovers (x2 at A1, x1 at A2).
    t1-injection of x_l' with letter beta: transport tau = pibar beta^-1 Dbar (tau M0 = anchor); stream frames tau(L0),
    tau(N1), tau(N2), tau(Y).  t2-injection with delta: tau' = pibar delta^-1 c2bar (tau' L0 = anchor)."""
    M0, N1, N2 = cf['M0'], cf['N1'], cf['N2']
    anchors = [cf['P1'], cf['Q']]
    out = []
    left = [(1, cf['A1']), (0, cf['A2'])]          # (array of the leftover, its global frame)
    for (lo, A) in left:
        pib = pi_bar(L, lo); pinv = sinv(pib)
        lam = dist(anchors[lo], A)
        X = image(pinv, A)                          # pibar^-1 A, at distance lam from L0
        mechs = []
        # t1 injections: tau(Y) = A for Y in {N2} U nbr(N2): beta(pibar^-1 A) = Dbar Y
        for Y in [N2] + cf['nN2']:
            if dist(M0, Y) == lam:
                mechs.append(('t1', 'N2' if Y == N2 else 'end', Y, image(DS, Y)))
        # t1 injections seen through t2: tau(Y') = A for Y' in {N1} U nbr(N1) (z on t2, or y on t2)
        for Y in [N1] + sorted(neighbors(N1)):
            if dist(M0, Y) == lam and dist(L0, Y) <= 2:
                mechs.append(('t1', 'N1' if Y == N1 else 'endt2', Y, image(DS, Y)))
        # t2 injections: tau'(Y) = A: delta(pibar^-1 A) = c2bar Y
        for Y in [N2] + cf['nN2']:
            if dist(L0, Y) == lam:
                mechs.append(('t2', 'N2' if Y == N2 else 'end', Y, image(L['c2'], Y)))
        rng.shuffle(mechs)
        if ONLY_N2:
            mechs = [m_ for m_ in mechs if m_[0] == 't1' and m_[1] == 'N2'] * n_per
        for (site, where, Y, Xp) in mechs[:max(1, n_per)]:
            got = 0; tries = 0
            while got < 1 and tries < 40:
                tries += 1
                b = find_P(rng, X, Xp, tries=200000)
                if b is None:
                    break
                if site == 't1':
                    tau = smul(pib, smul(sinv(b), DS))
                else:
                    tau = smul(pib, smul(sinv(b), L['c2']))
                o = sorder(tau)
                if odd_only and o % 2 == 0:
                    continue
                out.append(dict(site=site, l=lo, where=where, Y=Y, letter=b, tau_order=o, lam=lam))
                got += 1
    return out


# ------------------------------------------------------------------ program
class Prog:
    def __init__(self, E, S):
        self.E = E; self.S = S
        self.gens = {}            # name -> list of (l, gid, coefdict)

    def mulD(self, Mx, gid):
        return self.E.intern(self.E.snap(self.E.D @ Mx @ self.E.MATS[gid]))


def build_prog(E, Sy, L, injs):
    P = Prog(E, Sy)
    mats = {k: lift_P(E, v) for k, v in L.items()}
    for k, v in L.items():                                   # check the lifts
        gid = E.intern(E.snap(mats[k]))
        assert Sy.of(gid) == v, k
        assert E.is_mono(mats[k])
    ID = E.ID
    s = [P.mulD(mats['s1'], ID), P.mulD(mats['s2'], ID)]
    nal = 0
    t1 = [(0, P.mulD(mats['a1'], s[0]), {0: 1.0}), (1, P.mulD(mats['a2'], s[1]), {0: 1.0})]
    s2gen = None; P.copy_idx = None
    if 's2b' in L:                     # w7-b5: two-letter start s2 = D(sigma2 + lam sigma2b) x2, lam a free unknown
        nal += 1; P.copy_idx = nal
        s2b = P.mulD(mats['s2b'], ID)
        s2gen = [(1, s[1], {0: 1.0}), (1, s2b, {nal: 1.0})]
        t1.append((1, P.mulD(mats['a2'], s2b), {nal: 1.0}))
    if 'a2b' in L:                     # w7-b5: two-letter merge (a2 + lam a2b) s2
        nal += 1; P.copy_idx = nal
        t1.append((1, P.mulD(mats['a2b'], s[1]), {nal: 1.0}))
    t2i = []
    for j in injs:
        nal += 1; j['alpha'] = nal
        Bm = j['mat'] if 'mat' in j else lift_P(E, j['letter'])
        g = E.intern(E.snap(E.D @ Bm))
        if j['site'] == 't1':
            t1.append((j['l'], g, {nal: 1.0}))
        else:
            t2i.append((j['l'], g, {nal: 1.0}))
    t2 = [(l, P.mulD(mats['c2'], g), c) for (l, g, c) in t1] + t2i
    t3 = [(l, P.mulD(mats['c3'], g), c) for (l, g, c) in t2]
    P.gens = {'t1': merge(t1), 't2': merge(t2), 't3': merge(t3)}
    if s2gen is not None:
        P.gens['s2'] = merge(s2gen)
        s = [s[0], ID]                 # array 2's start is no longer a single term: only layer-0 / D b x2 rows are free
    P.s = s; P.nal = nal; P.mats = mats
    return P


def merge(T):
    d = {}
    for (l, g, c) in T:
        cc = d.setdefault((l, g), {})
        for k, v in c.items():
            cc[k] = cc.get(k, 0) + v
    return [(l, g, c) for (l, g), c in d.items()]


def clean_letter(P, L, l):
    """b0 in B with D b0 h pi_l in BK (the clean end letter for array l), lifted exactly."""
    E = P.E
    hpi = P.gens['t3'][[i for i, t in enumerate(P.gens['t3']) if t[0] == l and t[2] == {0: 1.0}][0]][1]
    X = image(P.S.of(hpi), Fstar)                     # (h pi)bar F*, a neighbour of L0
    rng = random.Random(7 + l)
    b = find_P(rng, X, image(DS, L0))
    Bm = lift_P(E, b)
    r = E.intern(E.snap(E.D @ Bm @ E.MATS[hpi]))
    assert P.S.frame(r) == Fstar
    return E.intern(E.snap(Bm))


def closure(P, target_l, dmin, b0, cap=3000):
    """component of the seed column y_{t3, b0} (exact, see NOTE): rows = non-free group elements of layer >= dmin;
    columns y_{u,b}, z_{u,b} (u in t1, t2, t3) touching them.  Returns rows, cols, entries, class columns (those with a
    term in B K x_target; only these carry kappa's target part)."""
    E, Sy = P.E, P.S
    Dinv = np.conj(E.D).T
    idx = {}
    for u, T in P.gens.items():
        for (l, g, c) in T:
            idx.setdefault((l, Sy.frame(g)), []).append((u, g))
    Fs = [Sy.frame(P.s[0]), Sy.frame(P.s[1])]
    FK = Fstar
    Dinv_cache = {}

    def dinv_frame(gid):
        v = Dinv_cache.get(gid)
        if v is None:
            h = E.intern(E.snap(Dinv @ E.MATS[gid])); v = (h, Sy.frame(h)); Dinv_cache[gid] = v
        return v

    def free(l, g):
        if E.layer(g) == 0:
            return True
        F = Sy.frame(g)
        if F == Fs[l] or F == FK:
            return True
        h, Fh = dinv_frame(g)
        if E.layer(h) == 0 or Fh == Fs[l]:
            return True
        return False

    rows = {}; queue = []; cols = {}; ent = []; klass = []; kterms = {}

    def add_col(typ, u, b):
        key = (typ, u, b)
        if key in cols:
            return True
        c = len(cols); cols[key] = c
        if c >= cap:
            return False
        bm = E.MATS[b]
        pre = E.D @ bm if typ == 'y' else bm
        iscl = False
        for (l2, g2, cf) in P.gens[u]:
            g3 = E.intern(E.snap(pre @ E.MATS[g2]))
            if l2 == target_l and E.layer(g3) == M and Sy.frame(g3) == FK:
                iscl = True
                kterms.setdefault(c, []).append(((l2, g3), cf))
            if E.layer(g3) < dmin:
                continue
            t3 = (l2, g3)
            r = rows.get(t3)
            if r is None:
                if free(l2, g3):
                    continue
                r = len(rows); rows[t3] = r; queue.append(t3)
            ent.append((r, c, cf))
        if iscl:
            klass.append(c)
        return True

    add_col('y', 't3', b0)
    while queue:
        (l, g) = queue.pop()
        cand = []
        for (u, gj) in idx.get((l, Sy.frame(g)), ()):
            cand.append(('z', u, E.MATS[g] @ np.conj(E.MATS[gj]).T))
        h, Fh = dinv_frame(g)
        for (u, gj) in idx.get((l, Fh), ()):
            cand.append(('y', u, E.MATS[h] @ np.conj(E.MATS[gj]).T))
        for (typ, u, bm) in cand:
            bm = E.snap(bm)
            if not E.is_mono(bm):
                continue
            if not add_col(typ, u, E.intern(bm)):
                return None
    return rows, cols, ent, (klass, kterms)


def assemble(cl, nal, rng):
    import scipy.sparse as sps
    rows, cols, ent, kl = cl
    klass, kterms = kl
    nr, nc = len(rows) + 1, len(cols)
    data = [([], [], []) for _ in range(nal + 1)]
    for (r, c, cf) in ent:
        for k, v in cf.items():
            data[k][0].append(r + 1); data[k][1].append(c); data[k][2].append(v)
    # normalisation row 0: a random functional of the K-PART (values on the B K x_target rows) equals 1
    vr = {}
    for c in klass:
        for (rk, cf) in kterms.get(c, []):
            if rk not in vr:
                vr[rk] = complex(rng.gauss(0, 1), rng.gauss(0, 1))
            for k, v in cf.items():
                data[k][0].append(0); data[k][1].append(c); data[k][2].append(vr[rk] * v)
    return [sps.csr_matrix((np.array(d[2], complex), (d[0], d[1])), shape=(nr, nc)) for d in data]


def vp_residual(Nss, alpha, want_J=True):
    rs, Js = [], []
    for Ns in Nss:
        A = Ns[0].toarray()
        for k in range(1, len(Ns)):
            if Ns[k].nnz:
                A += alpha[k - 1] * Ns[k].toarray()
        t = np.zeros(A.shape[0], complex); t[0] = 1
        if A.shape[1] == 0:
            rs.append(-t)
            if want_J:
                Js.append(np.zeros((A.shape[0], len(alpha)), complex))
            continue
        w, *_ = np.linalg.lstsq(A, t, rcond=1e-11)
        r = A @ w - t
        rs.append(r)
        if want_J:
            if len(Ns) > 1:
                V = np.stack([Ns[k] @ w for k in range(1, len(Ns))], axis=1)
                X = np.linalg.lstsq(A, V, rcond=1e-11)[0]
                Js.append(V - A @ X)
            else:
                Js.append(np.zeros((A.shape[0], 0), complex))
    r = np.concatenate(rs)
    if not want_J:
        return r, None
    return r, np.concatenate(Js, axis=0)


def lm(Nss, nal, rng, starts=3, iters=60):
    best = (np.inf, None)
    if nal == 0:
        r, _ = vp_residual(Nss, np.zeros(0, complex), False)
        return float(np.vdot(r, r).real), None
    for st in range(starts):
        al = np.array([rng.gauss(0, 1) + 1j * rng.gauss(0, 1) for _ in range(nal)])
        lamb = 1e-2
        r, J = vp_residual(Nss, al)
        cost = float(np.vdot(r, r).real)
        for it in range(iters):
            Jr = np.concatenate([np.concatenate([J.real, -J.imag], 1), np.concatenate([J.imag, J.real], 1)], 0)
            rr = np.concatenate([r.real, r.imag])
            g = Jr.T @ rr; H = Jr.T @ Jr
            improved = False
            for _ in range(10):
                try:
                    dx = -np.linalg.solve(H + lamb * (np.diag(np.diag(H)) + 1e-12 * np.eye(len(g))), g)
                except np.linalg.LinAlgError:
                    lamb *= 10; continue
                aln = al + dx[:nal] + 1j * dx[nal:]
                rn, _ = vp_residual(Nss, aln, False)
                cn = float(np.vdot(rn, rn).real)
                if cn < cost:
                    al, cost = aln, cn; lamb = max(lamb / 5, 1e-12); improved = True; break
                lamb *= 5
            if not improved or cost < 1e-24:
                break
            r, J = vp_residual(Nss, al)
        if cost < best[0]:
            best = (cost, al)
        if best[0] < 1e-20:
            break
    return best


def main():
    global ONLY_N2
    seed, lam1, lam2, nper, odd_only = [int(a) for a in sys.argv[1:6]]
    if odd_only == 2:
        ONLY_N2 = True; odd_only = 0
    dlist = [int(c) for c in sys.argv[6]] if len(sys.argv) > 6 else [5, 4, 3, 2, 1]
    rng = random.Random(seed); t0 = time.time()
    cf = make_config(rng, lam1, lam2)
    L = build_letters(rng, cf)
    injs = transports(rng, cf, L, nper, odd_only) if nper > 0 else []
    print("seed %d lam (%d,%d): %d injections offered: %s (%.0fs)" % (seed, lam1, lam2, len(injs), [
        (j['site'], 'x%d' % (j['l'] + 1), j['where'], 'lam%d' % j['lam'], 'ord%d' % j['tau_order']) for j in injs], time.time() - t0), flush=True)
    E = Eng(M, seed=seed); Sy = Symp(E)
    # sanity: symplectic model of layer and frames
    assert Sy.of(E.intern(E.D)) == DS and Sy.of(E.KID) == KS
    P = build_prog(E, Sy, L, injs)
    for u, T in P.gens.items():
        for (l, g, c) in T[:6]:
            assert E.layer(g) == dist(L0, Sy.frame(g))
    lay = {u: sorted(set((l, E.layer(g)) for (l, g, c) in T)) for u, T in P.gens.items()}
    print("  generator (array, layer) sets: %s" % lay, flush=True)
    b0s = [clean_letter(P, L, l) for l in (0, 1)]
    for d in dlist:
        cls = []
        for l in (0, 1):
            cl = closure(P, l, d, b0s[l])
            if cl is None:
                print("  d=%d: closure cap hit" % d, flush=True); return
            cls.append(cl)
        for li, cl in enumerate(cls):
            rows_, cols_, ent_, kl_ = cl
            prod = {}
            for (r, c, cf) in ent_:
                prod.setdefault(r, set()).add(c)
            inv_rows = {v: k for k, v in rows_.items()}
            single = [r for r, cs in prod.items() if len(cs) == 1]
            colkind = {}
            for k, c in cols_.items():
                colkind[c] = k[0] + '_' + k[1]
            byc = {}
            for r in single:
                c = next(iter(prod[r])); byc.setdefault(colkind[c], 0); byc[colkind[c]] += 1
            killed = set(next(iter(prod[r])) for r in single)
            print("     out %d: %d rows, %d cols; rows with ONE producer: %d (by producer kind %s); columns forced to 0 by them: %d of %d %s" % (
                li + 1, len(rows_), len(cols_), len(single), byc, len(killed), len(cols_),
                sorted(set(colkind[c] for c in killed))), flush=True)
        Nss = [assemble(cl, P.nal, rng) for cl in cls]
        cost, al = lm(Nss, P.nal, rng)
        print("  d=%d: closures (rows x cols, class cols) %s  residual %.3e  (%.0fs)" % (
            d, [(len(c[0]), len(c[1]), len(c[3][0])) for c in cls], cost, time.time() - t0), flush=True)
        if cost > 1e-12:
            r, _ = vp_residual(Nss, al if al is not None else np.zeros(0), False)
            n0 = Nss[0][0].shape[0]
            print("     per output: %.3e %.3e" % (float(np.vdot(r[:n0], r[:n0]).real), float(np.vdot(r[n0:], r[n0:]).real)), flush=True)
            break
    print("done (%.0fs)" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
