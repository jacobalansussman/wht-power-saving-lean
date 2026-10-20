"""b5eng.py (mw2-b5beyond): numerical LAYERED solver for W-in-M programs at block m (default 5), engine/widened letters.

Exact group elements from mw-big56/mwlib.Eng (2^m x 2^m dyadic matrices, interned).  A program is a list of moves;
move k has content sum_s alpha_s b_s u_{src_s} (src = earlier generator or an input x_l), u_k = D content.  Each
generator is a term list [(l, gid)] with a coefficient vector that is MULTILINEAR in alpha (forward pass).  The end
test for array l: [K]x_l in the C[B]-span of generators (Wd closure over aligned Z4-monomial images, mwlib), solved by
least squares (variable projection), RESTRICTED TO ROWS OF LAYER >= d (layer = log2 max row support = frame distance).
res_d = min_beta || (A beta - e_K)|_{layer >= d} ||^2 summed over arrays.  d = m: top layer only; d = 1: full unit.
LM (Levenberg-Marquardt) over alpha with the variable-projection Jacobian.
Labels: everything here is NUMERICAL evidence; a zero residual must be reconstructed and checked exactly (mwlib)."""
import sys, os, random, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'mw-big56'))
import numpy as np
from mwlib import Eng as _Eng0, SC, Z4
import hashlib


class Eng(_Eng0):
    """mwlib.Eng with compact storage: 16-byte digests as intern / frame keys, complex128 matrices (entries are
    multiples of 2^-14 of modulus <= 1)."""
    @staticmethod
    def key(M):
        return hashlib.blake2b(np.round(M.real * SC).astype(np.int32).tobytes() + np.round(M.imag * SC).astype(np.int32).tobytes(), digest_size=16).digest()

    def intern(self, M):
        k = self.key(M); i = self.IDX.get(k)
        if i is None:
            i = len(self.MATS); self.IDX[k] = i; self.MATS.append(M)
        return i

    def frame(self, gid):
        v = self._frm.get(gid)
        if v is None:
            M = self.MATS[gid]; rows = []
            for i in range(self.n):
                r = M[i]; nz = np.nonzero(np.abs(r) > 1e-9)[0]; c = r[nz[0]]
                for z in Z4:
                    a = np.angle(c * z)
                    if -1e-9 <= a < np.pi / 2 - 1e-9:
                        r = self.snap((r * z)[None, :])[0]; break
                rows.append(np.round(r.real * SC).astype(np.int32).tobytes() + np.round(r.imag * SC).astype(np.int32).tobytes())
            v = hashlib.blake2b(b"".join(sorted(rows)), digest_size=16).digest(); self._frm[gid] = v
        return v

    def closure_layered(self, gens, target, cap=4000, dmin=1):
        """mwlib closure restricted to rows of layer >= dmin (the projected identity on layers >= dmin needs no other
        rows: a column enters only through an aligned row of layer >= dmin)."""
        index = {}
        for i, g in enumerate(gens):
            for (l, gid) in g:
                if self.layer(gid) >= dmin:
                    index.setdefault((self.frame(gid), l), []).append((i, gid))
        rows = {target: 0}; cols = {}; ent = []; queue = [target]
        while queue:
            (lt, gt) = queue.pop()
            for (i, gterm) in index.get((self.frame(gt), lt), ()):
                bm = self.snap(self.MATS[gt] @ np.conj(self.MATS[gterm]).T)
                if not self.is_mono(bm):
                    continue
                b = self.intern(bm)
                if (i, b) in cols:
                    continue
                c = len(cols); cols[(i, b)] = c
                if c >= cap:
                    return None
                for (l2, g2), co in gens[i].items():
                    if self.layer(g2) < dmin:
                        continue
                    g3 = self.mulid(bm, g2)
                    t3 = (l2, g3)
                    r3 = rows.get(t3)
                    if r3 is None:
                        r3 = len(rows); rows[t3] = r3; queue.append(t3)
                    ent.append((r3, c, co))
        return rows, cols, ent

# ------------------------------------------------------------------ GL(m,2) helpers and letter menus
def gl_apply(P, x):
    y = 0; j = 0
    while x:
        if x & 1: y ^= P[j]
        x >>= 1; j += 1
    return y

def gl_mul(P, Q):
    return [gl_apply(P, q) for q in Q]

def gl_pow(P, k, m):
    R = [1 << j for j in range(m)]
    for _ in range(k): R = gl_mul(P, R)
    return R

def gl_order(P, m, cap=1000):
    I = [1 << j for j in range(m)]; R = list(P)
    for k in range(1, cap):
        if R == I: return k
        R = gl_mul(P, R)
    return None

def companion(poly_bits, d):
    """companion matrix of x^d + sum poly_bits[i] x^i on F_2^d: e_j -> e_{j+1}, e_{d-1} -> sum poly_bits[i] e_i"""
    P = [1 << (j + 1) for j in range(d - 1)]
    P.append(sum(1 << i for i in range(d) if poly_bits[i]))
    return P

def singer(m):
    polys = {2: [1, 1], 3: [1, 1, 0], 4: [1, 1, 0, 0], 5: [1, 0, 1, 0, 0], 6: [1, 1, 0, 0, 0, 0]}
    return companion(polys[m], m)

def block_diag(blocks, m):
    """blocks: list of (P_sub, d); acting on consecutive coordinate blocks"""
    P = []; off = 0
    for (Q, d) in blocks:
        for j in range(d): P.append(Q[j] << off)
        off += d
    while len(P) < m: P.append(1 << len(P))
    return P

def levi_menu(m):
    """odd-order Levi elements: Singer of F_2^m (order 2^m-1) and Singer cycles of sub-blocks (orders 3, 7, 15)"""
    out = {}
    eta = singer(m); out['eta%d' % gl_order(eta, m)] = eta
    for d in range(2, m):
        s = singer(d); o = gl_order(s, d)
        out['sub%d_%d' % (d, o)] = block_diag([(s, d)], m)
    if m >= 4:   # order 5 inside GL(4,2): Singer-15 cubed
        s4 = singer(4); out['sub4_5'] = block_diag([(gl_pow(s4, 3, 4), 4)], m)
    return out

class Letters:
    def __init__(self, E):
        self.E = E; m, n = E.m, E.n
        self.cache = {}

    def addr(self, P, c=0):
        return self.E.affine(A=list(P), c=c)

    def zp(self, b):
        n = self.E.n
        return np.diag([(-1) ** bin(b & x).count('1') for x in range(n)]).astype(complex)

    def sphase(self, a):
        """i^{a.x} (Clifford diagonal, linear part)"""
        n = self.E.n
        return np.diag([1j ** bin(a & x).count('1') for x in range(n)]).astype(complex)

    def cz(self, j, k):
        n = self.E.n
        return np.diag([(-1) ** (((x >> j) & (x >> k)) & 1) for x in range(n)]).astype(complex)


# ------------------------------------------------------------------ programs
class Prog:
    def __init__(self, E, W):
        self.E = E; self.W = W
        self.gens = [[(l, E.ID)] for l in range(W)]
        self.mvmaps = [None] * W        # per generator: list of (alpha index, src, si, ri) or None for inputs
        self.nal = 0; self.fixed = []
        self.desc = []

    def add_move(self, summ, Dm=None, desc=None):
        """summ: list of (src index, matrix b).  u = Dm . sum_s alpha_s b_s u_src.  Returns generator index."""
        E = self.E
        Dm = E.D if Dm is None else Dm
        terms = {}; order = []; maps = []
        for k, (src, b) in enumerate(summ):
            a = self.nal + k
            Db = E.snap(Dm @ b)
            si, ri = [], []
            for j, (l, g) in enumerate(self.gens[src]):
                g1 = E.intern(E.snap(Db @ E.MATS[g]))
                key = (l, g1)
                if key not in terms:
                    terms[key] = len(order); order.append(key)
                si.append(j); ri.append(terms[key])
            maps.append((a, src, np.array(si, dtype=np.int64), np.array(ri, dtype=np.int64)))
        self.fixed.append(self.nal)
        self.nal += len(summ)
        self.gens.append(order); self.mvmaps.append(maps); self.desc.append(desc)
        return len(self.gens) - 1

    def coefs(self, alpha):
        out = []
        for i, T in enumerate(self.gens):
            if self.mvmaps[i] is None:
                out.append(np.ones(1, dtype=complex)); continue
            v = np.zeros(len(T), dtype=complex)
            for (a, src, si, ri) in self.mvmaps[i]:
                np.add.at(v, ri, alpha[a] * out[src][si])
            out.append(v)
        return out

    def setup_final(self, cap=2000, maxcells=8e6, maxterms=2000, dmin=1):
        """memory guard: refuse programs whose dense closure matrices would exceed maxcells entries (16 B each)"""
        E = self.E
        if sum(len(T) for T in self.gens) > maxterms:
            self.why = 'terms'; return False
        tok = []; gd = []
        for i, T in enumerate(self.gens):
            d = {}
            for j, t in enumerate(T):
                d[t] = len(tok) + 1; tok.append((i, j))
            gd.append(d)
        sizes = [len(T) for T in self.gens]
        self.off = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
        self.fin = []
        for l in range(self.W):
            cl = E.closure_layered(gd, (l, E.KID), cap, dmin)
            if cl is None:
                self.why = 'cap'; return False
            rows, cols, ent = cl
            if len(rows) * max(1, len(cols)) > maxcells or len(rows) * 40 > maxcells:
                self.why = 'cells %d x %d' % (len(rows), len(cols)); return False
            lay = np.zeros(len(rows), dtype=np.int64)
            for t, r in rows.items():
                lay[r] = E.layer(t[1])
            Ea = np.array([(r, c, tok[co - 1][0], tok[co - 1][1]) for (r, c, co) in ent], dtype=np.int64).reshape(-1, 4)
            self.fin.append((len(rows), len(cols), Ea[:, 0], Ea[:, 1], self.off[Ea[:, 2]] + Ea[:, 3], lay))
        return True

    def flat(self, alpha):
        return np.concatenate(self.coefs(alpha))

    def vp(self, alpha, d, free_idx, want_J=True, eps=1e-7):
        fl = self.flat(alpha)
        if not np.isfinite(fl).all():
            return None, None
        Dv = []
        if want_J:
            for a in free_idx:
                al = alpha.copy(); al[a] += eps
                Dv.append((self.flat(al) - fl) / eps)
        rs, Js = [], []
        for (nr, nc, R, C, F, lay) in self.fin:
            keep = lay >= d
            idx = -np.ones(nr, dtype=np.int64); idx[keep] = np.arange(keep.sum())
            sel = keep[R]
            R2, C2, F2 = idx[R[sel]], C[sel], F[sel]
            nr2 = int(keep.sum())
            A = np.zeros((nr2, nc), dtype=complex)
            np.add.at(A, (R2, C2), fl[F2])
            t = np.zeros(nr2, dtype=complex); t[idx[0]] = 1
            if nc == 0:
                rs.append(-t)
                if want_J: Js.append(np.zeros((nr2, len(Dv)), dtype=complex))
                continue
            beta = np.linalg.lstsq(A, t, rcond=1e-10)[0]
            r = A @ beta - t
            rs.append(r)
            if want_J:
                V = np.zeros((nr2, len(Dv)), dtype=complex)
                for k, dd in enumerate(Dv):
                    w = dd[F2] * beta[C2]
                    V[:, k] = np.bincount(R2, weights=w.real, minlength=nr2) + 1j * np.bincount(R2, weights=w.imag, minlength=nr2)
                X = np.linalg.lstsq(A, V, rcond=1e-10)[0]
                Js.append(V - A @ X)
        r = np.concatenate(rs)
        if not want_J:
            return r, None
        Jc = np.concatenate(Js, axis=0)
        Jr = np.concatenate([np.concatenate([Jc.real, -Jc.imag], axis=1),
                             np.concatenate([Jc.imag, Jc.real], axis=1)], axis=0)
        return r, Jr

    def lm(self, d, x0, free_idx, iters=60, tol=1e-24):
        N = self.nal; nf = len(free_idx)
        def full(x):
            al = np.ones(N, dtype=complex); al[free_idx] = x[:nf] + 1j * x[nf:]; return al
        x = x0.copy(); lam = 1e-3
        r, J = self.vp(full(x), d, free_idx)
        if r is None:
            return x, np.inf
        rr = np.concatenate([r.real, r.imag]); cost = rr @ rr
        for it in range(iters):
            g = J.T @ rr; Hm = J.T @ J
            improved = False
            for _ in range(12):
                try:
                    dx = -np.linalg.solve(Hm + lam * (np.diag(np.diag(Hm)) + 1e-12 * np.eye(len(x))), g)
                except np.linalg.LinAlgError:
                    lam *= 10; continue
                xn = x + dx
                rn, _ = self.vp(full(xn), d, free_idx, want_J=False)
                if rn is None:
                    lam *= 10; continue
                rrn = np.concatenate([rn.real, rn.imag]); cn = rrn @ rrn
                if cn < cost:
                    x, cost = xn, cn; lam = max(lam / 5, 1e-12); improved = True; break
                lam *= 5
            if not improved or cost < tol:
                break
            r, J = self.vp(full(x), d, free_idx)
            rr = np.concatenate([r.real, r.imag])
        return x, cost

    def free_idx(self):
        fx = set(self.fixed)
        return [a for a in range(self.nal) if a not in fx]

    def full_alpha(self, x):
        fi = self.free_idx(); nf = len(fi)
        al = np.ones(self.nal, dtype=complex); al[fi] = x[:nf] + 1j * x[nf:]
        return al


def plain_chain_test(m=5, W=2):
    """control: W separate plain chains (directions e1..em via coordinate transpositions) must give res 0 at every d."""
    E = Eng(m, seed=1); Lt = Letters(E)
    P = Prog(E, W)
    for l in range(W):
        src = l
        for j in range(m):
            # u = D_{e_{j+1}} u_prev realised as D conjugated: dir(e_{j+1}) = S D S^-1 with S swapping e1 <-> e_{j+1}
            if j == 0:
                b = E.I
                src = P.add_move([(src, b)])
            else:
                S = [1 << t for t in range(m)]; S[0], S[j] = S[j], S[0]
                Sm = Lt.addr(S)
                # D S u_prev = S dir(e_{j+1}) u_prev; next content multiplies by S^-1 = S: keep S on the left
                src = P.add_move([(src, Sm)])
                # undo by an extra factor at the next move (S is an involution): fold into readout via closure
    ok = P.setup_final()
    return P, ok
