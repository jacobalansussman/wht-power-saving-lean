"""s4lib.py (mw2-singer4): exact register-form words at block m (default 4) with the Singer element of F_2^m/<x>,
and a JOINT end solver: W outputs, W shared top moves E_j c_j with UNKNOWN contents c_j (solved), fixed separation
coefficients alpha_{o,j} in C[<eta>] (step 1 classes), plus the prefix span; exact (closure over frames, rank mod p in C,
then an exact rational solve only if mod p passes).  Reuses tb4-beam bm.Game (exact 2^m x 2^m dyadic matrices) and
mw-sym4 symlib.Sym.hmax_elem (Levi lifts c_P with K^-1 c_P^-1 K free).

Generators: dict {(l, gid): Fraction}.  Free group B' = every monomial with entries in {1, i, -1, -i} (widened game;
the engine letters used to BUILD prefixes are affine maps, Clifford phases, Paulis and the Levi lifts, all in E).
Row of the joint system: (o, l, gid) = coefficient of [g] x_l in output o; target (o, o, K) = 1, all else 0; rows at free
g are dropped (raw inputs absorb them).  Columns: ('S', o, i, b) = output o uses b * gen_i;  ('E', j, i, b) = content
of top j contains b * gen_i (shared by every output o with alpha_{o,j} != 0)."""
import sys, os, os, itertools, random, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../common'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../symmetric-block4-skeletons'))
import numpy as np
from fractions import Fraction as Fr
import bm as BM
from bm import Game, Z4
from symlib import Sym, gl_apply

SC = 4096.0

def fast_frame_key(M):
    """key of the left coset B'M (B' = Z4-monomials): normalise each row's phase by its first nonzero entry
    (angle in [0, pi/2)), round, sort rows."""
    nz = np.abs(M) > 1e-9
    first = np.argmax(nz, axis=1)
    v = M[np.arange(M.shape[0]), first]
    ang = np.angle(v)
    k = np.mod(-np.floor(ang / (np.pi / 2) + 1e-9), 4).astype(int)
    rot = np.array([1, 1j, -1, -1j])[k]
    R = M * rot[:, None]
    I = np.concatenate([np.round(R.real * SC), np.round(R.imag * SC)], axis=1).astype(np.int64)
    order = np.lexsort(I.T[::-1])
    return I[order].tobytes()

class S4:
    def __init__(self, m=4, seed=1):
        self.m = m; self.n = 1 << m
        self.sym = Sym(m); self.G = G = self.sym.G
        self.rng = random.Random(seed)
        self.FK = {}; self.MONO = {}; self.INV = {}
        self.D1 = G.D1
        # Singer element of F_2^m/<e1> (m = 4: companion of t^3 + t + 1 on e2, e3, e4), fixing e1
        if m == 4 and os.environ.get('S4_ORDER', '7') == '3':
            P = [1, 4, 8, 2]               # order 3: cyclic shift e2 -> e3 -> e4 -> e2 (v-ms-fano's sigma), fixing e1
        elif m == 4:
            P = [1, 4, 8, 2 | 4]
        elif m == 3:
            P = [1, 4, 2 | 4]              # order 3 on F_2^2
        else:
            raise ValueError
        self.P = P
        c, cp = self.sym.hmax_elem(P)
        self.eta = c; self.eta_id = G.intern(c)
        o = 1; M = c
        while not np.allclose(M, np.eye(self.n)) and o < 100:
            M = G.mul(M, c); o += 1
        self.eta_order_proj = o
        # powers of eta (exact)
        self.etap = [np.eye(self.n, dtype=complex)]
        for _ in range(1, o if o in (3, 7) else (7 if m == 4 else 3)):
            self.etap.append(G.mul(self.etap[-1], c))
        self.etap_id = [G.intern(x) for x in self.etap]

    # ---------------------------------------------------------------- cached group facts
    def fk(self, gid):
        v = self.FK.get(gid)
        if v is None:
            v = fast_frame_key(self.G.MATS[gid]); self.FK[gid] = v
        return v

    def inv(self, gid):
        v = self.INV.get(gid)
        if v is None:
            v = self.G.intern(np.conj(self.G.MATS[gid]).T); self.INV[gid] = v
        return v

    def mulid(self, *ids):
        G = self.G
        return G.intern(G.mul(*[G.MATS[i] for i in ids]))

    def prod(self, *ids):
        """product matrix WITHOUT interning (memory: only monomials and kept terms are interned)"""
        G = self.G
        return G.mul(*[G.MATS[i] for i in ids])

    def mono_id(self, *ids):
        """intern the product only if it is a monomial; else None"""
        M = self.prod(*ids)
        if not self.G.is_monomial(M):
            return None
        return self.G.intern(M)

    def is_mono(self, gid):
        v = self.MONO.get(gid)
        if v is None:
            v = self.G.is_monomial(self.G.MATS[gid]); self.MONO[gid] = v
        return v

    def depth(self, gid):
        return int(round(np.log2(self.G.rowsupp(self.G.MATS[gid]).max())))

    def dK(self, gid):
        return self.G.dist(gid)

    # ---------------------------------------------------------------- engine free letters
    def rand_free(self, pauli=True):
        m, n = self.m, self.n; rng = self.rng
        while True:
            A = [rng.randrange(1, n) for _ in range(m)]
            span = {0}; ok = True
            for a in A:
                if a in span:
                    ok = False; break
                span |= {s ^ a for s in span}
            if ok:
                break
        c = rng.randrange(n)
        sig = []
        for x in range(n):
            y = c
            for j in range(m):
                if (x >> j) & 1:
                    y ^= A[j]
            sig.append(y)
        a = [rng.randrange(4) for _ in range(m)]
        B = {(j, k): rng.randrange(2) for j in range(m) for k in range(j + 1, m)}
        q = []
        for x in range(n):
            b = [(x >> j) & 1 for j in range(m)]
            v = sum(a[j] * b[j] for j in range(m)) + 2 * sum(B[(j, k)] * b[j] * b[k] for j in range(m) for k in range(j + 1, m))
            if pauli:
                z = 0
            q.append(v % 4)
        if pauli:
            z = rng.randrange(n)
            q = [(q[x] + 2 * (bin(z & x).count('1') & 1)) % 4 for x in range(n)]
        return self.G.mono(sig, q)

    # ---------------------------------------------------------------- generators
    def inp(self, l):
        return {(l, self.G.ID): Fr(1)}

    def lmul(self, M, gen):
        G = self.G; out = {}
        for (l, g), c in gen.items():
            t = (l, G.intern(G.mul(M, G.MATS[g])))
            out[t] = out.get(t, 0) + c
        return {t: c for t, c in out.items() if c != 0}

    @staticmethod
    def add(*pairs):
        out = {}
        for (c0, gen) in pairs:
            for t, c in gen.items():
                out[t] = out.get(t, 0) + c0 * c
        return {t: c for t, c in out.items() if c != 0}

    def move(self, summands, Dm=None):
        """summands: list of (coef, monomial matrix or None, gen); returns D * sum"""
        G = self.G
        parts = []
        for (co, b, gen) in summands:
            parts.append((Fr(co), gen if b is None else self.lmul(b, gen)))
        c = self.add(*parts)
        return self.lmul(self.D1 if Dm is None else Dm, c)

    # ---------------------------------------------------------------- joint end solver
    def joint(self, gens, W, Es, alpha, cap=2500, rowcap=12000, exact=True):
        """gens: prefix generators (inputs need not be listed).  Es: list of W paid matrices (rank-1 moves D1 b).
        alpha[o][j] = dict {k: coef} meaning sum_k coef * eta^k.  Returns ('yes', info) / ('no', stats) / ('cap', stats)."""
        G = self.G
        index = {}
        for i, g in enumerate(gens):
            for (l, gid) in g:
                index.setdefault((self.fk(gid), l), []).append((i, gid))
        Eids = [G.intern(E) for E in Es]
        Einv = [self.inv(e) for e in Eids]
        etainv = [self.inv(e) for e in self.etap_id]
        # pre-multipliers for top j, power k: eta^k E_j
        pre = {}
        for j in range(W):
            for o in range(W):
                for k in alpha[o][j]:
                    if (j, k) not in pre:
                        pre[(j, k)] = self.mulid(self.etap_id[k], Eids[j])
        rows = {}; cols = {}; ent = []; queue = []
        def addrow(t):
            r = rows.get(t)
            if r is None:
                if len(rows) >= rowcap:
                    raise OverflowError
                r = len(rows); rows[t] = r; queue.append(t)
            return r
        try:
            for o in range(W):
                addrow((o, o, G.KID))
            while queue:
                (o, l, g) = queue.pop()
                cands = []
                for (i, gt) in index.get((self.fk(g), l), ()):
                    cands.append((('S', o, i), self.mulid(g, self.inv(gt))))
                for j in range(W):
                    for k in alpha[o][j]:
                        h = self.mulid(Einv[j], etainv[k], g)       # E_j^-1 eta^-k g
                        for (i, gt) in index.get((self.fk(h), l), ()):
                            cands.append((('E', j, i), self.mulid(h, self.inv(gt))))
                for (key, b) in cands:
                    if not self.is_mono(b):
                        continue
                    ck = key + (b,)
                    if ck in cols:
                        continue
                    c = len(cols); cols[ck] = c
                    if c >= cap:
                        return 'cap', {'rows': len(rows), 'cols': len(cols)}
                    if key[0] == 'S':
                        _, oo, i = key
                        for (l2, g2), co in gens[i].items():
                            g3 = self.mulid(b, g2)
                            if self.is_mono(g3):
                                continue
                            ent.append((addrow((oo, l2, g3)), c, co))
                    else:
                        _, j, i = key
                        for o2 in range(W):
                            for k, ak in alpha[o2][j].items():
                                for (l2, g2), co in gens[i].items():
                                    g3 = self.mulid(pre[(j, k)], b, g2)
                                    if self.is_mono(g3):
                                        continue
                                    ent.append((addrow((o2, l2, g3)), c, co * ak))
        except OverflowError:
            return 'cap', {'rows': len(rows), 'cols': len(cols)}
        nr, nc = len(rows), len(cols)
        if nc == 0:
            return 'no', {'rows': nr, 'cols': 0}
        P = BM.P
        A = np.zeros((nr, nc + 1))
        for (r, c, co) in ent:
            co = Fr(co)
            A[r, c] = (A[r, c] + (co.numerator % P) * pow(co.denominator % P, P - 2, P)) % P
        for o in range(W):
            A[rows[(o, o, G.KID)], nc] = 1
        Rr, piv = BM.rref(A, full=False)
        st = {'rows': nr, 'cols': nc, 'rank': len(piv)}
        if nc in piv:
            return 'no', st
        if not exact:
            return 'yes-modp', st
        x = self.solve_Q(nr, nc, ent, [rows[(o, o, G.KID)] for o in range(W)])
        if x is None:
            return 'no-Q', st
        inv = {c: k for k, c in cols.items()}
        return 'yes', {'st': st, 'sol': [(inv[c], x[c]) for c in range(nc) if x[c]]}

    @staticmethod
    def solve_Q(nr, nc, ent, trows):
        """sparse exact elimination over Q (dict rows)"""
        R = [dict() for _ in range(nr)]
        for (r, c, co) in ent:
            R[r][c] = R[r].get(c, Fr(0)) + Fr(co)
        rhs = [Fr(0)] * nr
        for t in trows:
            rhs[t] = Fr(1)
        rowsl = [(R[r], rhs[r]) for r in range(nr)]
        pivots = {}
        done = []
        for (row, b) in rowsl:
            row = {c: v for c, v in row.items() if v != 0}
            # reduce by existing pivots
            changed = True
            while True:
                cs = [c for c in row if c in pivots]
                if not cs:
                    break
                c = cs[0]
                prow, pb = pivots[c]
                f = row[c]
                for cc, vv in prow.items():
                    nv = row.get(cc, Fr(0)) - f * vv
                    if nv == 0:
                        row.pop(cc, None)
                    else:
                        row[cc] = nv
                b = b - f * pb
            if not row:
                if b != 0:
                    return None
                continue
            c = min(row)
            f = row[c]
            row = {cc: vv / f for cc, vv in row.items()}
            b = b / f
            # eliminate c from other pivots (keep reduced form lazily: not needed for consistency)
            pivots[c] = (row, b)
        # back substitution: free vars 0
        x = [Fr(0)] * nc
        order = sorted(pivots, reverse=True)
        for c in order:
            row, b = pivots[c]
            s = b
            for cc, vv in row.items():
                if cc != c:
                    s -= vv * x[cc]
            x[c] = s
        return x

    # ---------------------------------------------------------------- exact verification of a finished word
    def verify(self, gens, W, outputs):
        """outputs[o] = list of (coef, monomial gid, generator index or ('in', l)); checks sum == K x_o exactly"""
        G = self.G
        for o in range(W):
            acc = {}
            for (co, b, i) in outputs[o]:
                gen = self.inp(i[1]) if isinstance(i, tuple) else gens[i]
                for (l, g), c in gen.items():
                    t = (l, self.mulid(b, g))
                    acc[t] = acc.get(t, 0) + co * c
            acc = {t: c for t, c in acc.items() if c != 0}
            if acc != {(o, G.KID): 1}:
                return False, o, acc
        return True, None, None

def E_pool(S, gens, maxn=12, rng=None):
    """paid rank-1 elements E = K g^-1 for prefix terms g one move from the kernel frame, legal (E = b D1 b', bm.one_move)"""
    G = S.G; out = {}
    for gen in gens:
        for (l, gid) in gen:
            if G.dist(gid) != 1:
                continue
            E = G.mul(G.K, np.conj(G.MATS[gid]).T)
            k = fast_frame_key(E.T.copy())          # E up to a RIGHT free factor gives the same move set
            if k in out:
                continue
            if not (G.rowsupp(E) == 2).all() or not G.one_move(E):
                continue
            out[k] = E
    vals = list(out.values())
    if rng is not None:
        rng.shuffle(vals)
    return vals[:maxn]
