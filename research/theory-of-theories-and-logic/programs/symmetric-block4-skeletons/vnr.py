"""vnr.py -- v-nag-rerun (skeptic): MY OWN rebuild of a symmetric skeleton from its dumped letters (recipes/NAME.json),
written from the semantics stated in mw-sym4 symsolve.add_move / bm.Game.closure (read, not imported), with EVERY
coefficient free (no normalisation, no charts, no forced zeros).

Exact arithmetic: matrices are int64 arrays scaled by S = 4096 (products divided back exactly, asserted); frame key =
coset B g (B = Z4-monomials on the left): each row rotated by the unique i^k making its first nonzero entry have
re > 0, im >= 0, rows then sorted (my own canonical form, different from v-nag's min-over-rotations).
Polynomials: dict sorted-tuple-of-alpha-indices -> int (characters here are +-1; asserted).

Outputs (functions): relaxation rows, full closure system for an array, Singular exports, mod-p screens."""
import json, sys, random
import numpy as np
from fractions import Fraction

S = 4096

class Mats:
    def __init__(self, m):
        self.m = m; n = self.n = 1 << m
        I = np.eye(n, dtype=np.int64) * S; Z = np.zeros((n, n), dtype=np.int64)
        self.I = (I, Z)
        self.D = {}
        for z in range(1, n):
            X = np.zeros((n, n), dtype=np.int64)
            for x in range(n):
                X[x ^ z, x] = S
            self.D[z] = ((I + X) // 2, (I - X) // 2)       # (1+i)/2 I + (1-i)/2 X
        K = self.I
        for j in range(m):
            K = self.mul(K, self.D[1 << j])
        self.K = K
        self.store = {}

    def mul(self, A, B):
        r = A[0] @ B[0] - A[1] @ B[1]; i = A[0] @ B[1] + A[1] @ B[0]
        assert not (r % S).any() and not (i % S).any()
        return (r // S, i // S)

    @staticmethod
    def inv(A):        # unitary: conjugate transpose
        return (A[0].T.copy(), -A[1].T.copy())

    def key(self, A):
        k = A[0].tobytes() + A[1].tobytes()
        if k not in self.store:
            self.store[k] = A
        return k

    def is_mono(self, A):
        nz = (A[0] != 0) | (A[1] != 0)
        if not ((nz.sum(0) == 1).all() and (nz.sum(1) == 1).all()):
            return False
        for a, b in zip(A[0][nz], A[1][nz]):
            if not ((abs(a) == S and b == 0) or (a == 0 and abs(b) == S)):
                return False
        return True

    def frame(self, A):
        R, Im = A; rows = []
        for x in range(self.n):
            r, i = R[x].tolist(), Im[x].tolist()
            j = next(k for k in range(self.n) if r[k] or i[k])
            # rotate by i^t until entry j has re > 0, im >= 0
            for _ in range(4):
                if r[j] > 0 and i[j] >= 0:
                    break
                r, i = [-v for v in i], r[:]     # multiply by i: (r + i im) * i = -im + i r
            rows.append(tuple(r) + tuple(i))
        rows.sort()
        return tuple(rows)

def lm(M, e):
    return (np.array(e[0], dtype=np.int64) * (S // 256), np.array(e[1], dtype=np.int64) * (S // 256))

def padd(P, Q, c=1):
    for k, v in Q.items():
        w = P.get(k, 0) + c * v
        if w == 0: P.pop(k, None)
        else: P[k] = w
    return P

def ptimes(P, a, c):
    out = {}
    for k, v in P.items():
        k2 = tuple(sorted(k + (a,))) if a is not None else k
        out[k2] = out.get(k2, 0) + c * v
    return {k: v for k, v in out.items() if v}

class Skel:
    def __init__(self, rec, zero=(), one=(), keptj=None):
        """keptj: per move, the list of tau indices j whose images are generators (from the all-free build); when given,
        exactly those images are added (no dedupe), so generator indices (the recipes' src) stay aligned."""
        """zero: alpha indices set to 0 (summand dropped).  Every other alpha is a free variable."""
        self.rec = rec; m = rec['m']; W = self.W = rec['W']
        M = self.M = Mats(m)
        self.H = []
        for (c, pi) in rec['Hel']:
            c = lm(M, c); assert M.is_mono(c)
            cp = M.mul(M.mul(M.inv(M.K), M.inv(c)), M.K)     # c' = K^-1 c^-1 K
            assert M.is_mono(cp)
            self.H.append((c, cp, tuple(pi)))
        nH = len(self.H)
        self.gens = []
        Ik = M.key(M.I)
        for l in range(W):
            self.gens.append({(l, Ik): {(): 1}})
        self.na = 0; self.cost = 0; self.gmove = [None] * W
        for t, (kind, summ, Dm) in enumerate(rec['moves']):
            Dmat = M.D[1] if Dm is None else lm(M, Dm)
            osum = kind in ('S', 'A', 'B', 'C', 'N', 'L')
            r = {'S': 0, 'A': 1, 'B': 2, 'C': 3}.get(kind, 0)
            Hs = list(range(nH)) if osum else [0]
            u = {}
            for k, (src, b) in enumerate(summ):
                a = self.na + k
                if a in zero:
                    continue
                avar = None if a in one else a
                b = lm(M, b); assert M.is_mono(b)
                for j in Hs:
                    ch = np.exp(2j * np.pi * r * j / nH)
                    chr_ = int(round(ch.real)); assert abs(ch - chr_) < 1e-12, "non-real character"
                    c, cp, pi = self.H[j]
                    for (l, gk), P in self.gens[src].items():
                        g1 = M.mul(b, M.store[gk])
                        if j == 0: l2, g2 = l, g1
                        else: l2, g2 = pi[l], M.mul(M.mul(c, g1), cp)
                        g3 = M.mul(Dmat, g2)
                        kk = (l2, M.key(g3))
                        u.setdefault(kk, {})
                        padd(u[kk], ptimes(P, avar, chr_))
            u = {k: P for k, P in u.items() if P}
            self.na += len(summ)
            reps = [u]; self.gens.append(u); self.gmove.append(t)
            nimg = 0
            self.keptj = getattr(self, 'keptj', [])
            kj = []
            for j in range(1, nH):
                img = self.tau(j, u)
                if keptj is not None:
                    if j not in keptj[t]: continue
                elif any(self.prop(img, w) for w in reps):
                    continue
                kj.append(j)
                reps.append(img); self.gens.append(img); self.gmove.append(t); nimg += 1
            self.keptj.append(kj)
            self.cost += 1 + nimg

    def tau(self, j, u):
        M = self.M; c, cp, pi = self.H[j]
        return {(pi[l], M.key(M.mul(M.mul(c, M.store[gk]), cp))): dict(P) for (l, gk), P in u.items()}

    def prop(self, v, w):
        """v == (monomial b) . w . (nonzero constant), exactly"""
        M = self.M
        if len(v) != len(w): return False
        if not w: return True
        (l0, g0), P0 = next(iter(w.items()))
        f0 = M.frame(M.store[g0])
        for (l1, g1), Q1 in v.items():
            if l1 != l0 or M.frame(M.store[g1]) != f0: continue
            b = M.mul(M.store[g1], M.inv(M.store[g0]))
            if not M.is_mono(b): continue
            ratio = None; ok = True
            for (l, gk), P in w.items():
                Q = v.get((l, M.key(M.mul(b, M.store[gk]))))
                if Q is None or set(Q) != set(P): ok = False; break
                for mo in P:
                    if ratio is None: ratio = Fraction(Q[mo], P[mo])
                    elif Fraction(Q[mo], P[mo]) != ratio: ok = False; break
                if not ok: break
            if ok: return True
        return False

    # ---------------------------------------------------------------- relaxation (Lemma 1)
    def relaxation(self, l0=0):
        M = self.M; rows = {}; fc = {}
        for s, u in enumerate(self.gens):
            for (l, gk), P in u.items():
                f = fc.get(gk)
                if f is None: f = fc[gk] = M.frame(M.store[gk])
                R = rows.setdefault((f, l), {})
                R[s] = padd(R.get(s, {}), P)
        out = []
        Kf = (M.frame(M.K), l0)
        for key, R in rows.items():
            R = {s: P for s, P in R.items() if P}
            const = -1 if key == Kf else 0
            if R or const: out.append((R, const))
        if Kf not in rows: out.append(({}, -1))
        return out

    # ---------------------------------------------------------------- full system (closure, f_s in C[B])
    def closure(self, l0=0, cap=200000):
        M = self.M
        fc = {}
        def fr(gk):
            f = fc.get(gk)
            if f is None: f = fc[gk] = M.frame(M.store[gk])
            return f
        index = {}
        for i, u in enumerate(self.gens):
            for (l, gk) in u:
                index.setdefault((fr(gk), l), []).append((i, gk))
        free = {}
        def isfree(gk):
            v = free.get(gk)
            if v is None: v = free[gk] = M.is_mono(M.store[gk])
            return v
        Kk = M.key(M.K)
        rows = {(l0, Kk): 0}; cols = {}; ent = {}; queue = [(l0, Kk)]
        while queue:
            lt, gt = queue.pop()
            for (i, gterm) in index.get((fr(gt), lt), ()):
                bm = M.mul(M.store[gt], M.inv(M.store[gterm]))
                assert M.is_mono(bm)
                bk = M.key(bm)
                if (i, bk) in cols: continue
                c = cols[(i, bk)] = len(cols)
                if c >= cap: return None
                for (l2, g2), P in self.gens[i].items():
                    g3 = M.key(M.mul(bm, M.store[g2]))
                    if isfree(g3): continue
                    t3 = (l2, g3)
                    r3 = rows.get(t3)
                    if r3 is None:
                        r3 = rows[t3] = len(rows); queue.append(t3)
                    E = ent.setdefault((r3, c), {})
                    padd(E, P)
        A = [dict() for _ in range(len(rows))]
        for (r, c), P in ent.items():
            if P: A[r][c] = P
        return A, len(cols)

# ---------------------------------------------------------------- polynomial helpers (Fraction coefficients)
def pmul(P, Q):
    out = {}
    for a, x in P.items():
        for b, y in Q.items():
            k = tuple(sorted(a + b)); w = out.get(k, 0) + x * y
            if w == 0: out.pop(k, None)
            else: out[k] = w
    return out

def const_eliminate(A, ncols, maxsteps=10**7):
    """A: list of rows {col: poly}; row 0 = target (rhs 1).  Represent rhs as column 'r'.
    Pivot only on entries that are nonzero CONSTANTS (valid for every alpha).  Returns the core rows."""
    rows = []
    for i, R in enumerate(A):
        R2 = {c: {k: Fraction(v) for k, v in P.items()} for c, P in R.items()}
        if i == 0: R2['r'] = {(): Fraction(-1)}
        if R2: rows.append(R2)
    alive = set(range(len(rows)))
    colrows = {}
    for i, R in enumerate(rows):
        for c in R:
            if c != 'r': colrows.setdefault(c, set()).add(i)
    while True:
        piv = None
        best = None
        for i in alive:
            R = rows[i]
            for c, P in R.items():
                if c != 'r' and len(P) == 1 and () in P:
                    sc = len(colrows.get(c, ())) * len(R)
                    if best is None or sc < best:
                        best = sc; piv = (i, c)
            if best is not None and best <= 2:
                break
        if piv is None: break
        i, c = piv
        Ri = rows[i]; ci = Ri[c][()]
        alive.discard(i)
        for cc in Ri:
            if cc != 'r': colrows[cc].discard(i)
        for j in list(colrows.get(c, ())):
            if j == i: continue
            R = rows[j]
            f = {k: -v / ci for k, v in R[c].items()}
            for cc, P in Ri.items():
                newP = dict(R.get(cc, {}))
                for k, v in pmul(f, P).items():
                    w = newP.get(k, 0) + v
                    if w == 0: newP.pop(k, None)
                    else: newP[k] = w
                if newP:
                    if cc not in R and cc != 'r': colrows.setdefault(cc, set()).add(j)
                    R[cc] = newP
                else:
                    if cc in R:
                        del R[cc]
                        if cc != 'r': colrows[cc].discard(j)
            assert c not in R
            if not R: alive.discard(j)
        colrows.pop(c, None)
    core = [rows[i] for i in sorted(alive) if rows[i]]
    return core

def sing_poly(P, aname, extra=None):
    terms = []
    for k, v in P.items():
        mon = []
        cnt = {}
        for a in k: cnt[a] = cnt.get(a, 0) + 1
        for a, e in sorted(cnt.items()):
            mon.append(aname(a) + (f"^{e}" if e > 1 else ""))
        if extra: mon.append(extra)
        v = Fraction(v)
        cs = f"({v.numerator}/{v.denominator})" if v.denominator != 1 else f"({v.numerator})"
        terms.append(cs + ("*" + "*".join(mon) if mon else ""))
    return "+".join(terms) if terms else "0"

def singular_text(rows_polys, avars, bvars, extra_ideal=(), tag='', char=0):
    """rows_polys: list of strings (already Singular syntax)."""
    allv = [f"a{a}" for a in avars] + [f"b{b}" for b in bvars]
    if not allv: allv = ['zz']
    s = f"ring r = {char},({','.join(allv)}),dp;\noption(redSB);\nideal I =\n" + ",\n".join(list(rows_polys) + list(extra_ideal)) + ";\n"
    s += "ideal G = std(I);\n"
    s += f'if (size(G) == 1 && G[1] == 1) {{ print("RESULT {tag} GB=1"); }} else {{ print("RESULT {tag} GB!=1 size " + string(size(G)) + " dim " + string(dim(G))); }}\nquit;\n'
    return s

def relax_singular(sk, tag, extra=()):
    rel = sk.relaxation()
    polys = []; av = set(); bv = set()
    for R, const in rel:
        P = {}
        for s, Q in R.items():
            for k, v in Q.items():
                P[k + (('b', s),)] = P.get(k + (('b', s),), 0) + v
                av.update(k); bv.add(s)
        terms = []
        for k, v in P.items():
            if v == 0: continue
            al = [x for x in k if not isinstance(x, tuple)]; bb = [x for x in k if isinstance(x, tuple)]
            cnt = {}
            for a in al: cnt[a] = cnt.get(a, 0) + 1
            mon = [f"a{a}" + (f"^{e}" if e > 1 else "") for a, e in sorted(cnt.items())] + [f"b{x[1]}" for x in bb]
            terms.append(f"({v})*" + "*".join(mon))
        if const: terms.append(f"({const})")
        polys.append("+".join(terms) if terms else "0")
    return singular_text(polys, sorted(av), sorted(bv), extra, tag), sorted(av), sorted(bv)

def core_singular(core, tag, extra=()):
    av = sorted({a for R in core for P in R.values() for k in P for a in k})
    bv = sorted({c for R in core for c in R if c != 'r'})
    polys = []
    for R in core:
        parts = []
        for c, P in R.items():
            parts.append(sing_poly(P, lambda a: f"a{a}", None if c == 'r' else f"b{c}"))
        polys.append("+".join(parts))
    return singular_text(polys, av, bv, extra, tag), av, bv

# ---------------------------------------------------------------- mod-p screens
PR = 1073741789     # = 1 mod 4 ?  (checked at runtime)
def rank_mod(A, p):
    A = A.copy() % p; r = 0; R, C = A.shape
    for c in range(C):
        if r == R: break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        pr = r + nz[0]
        if pr != r: A[[r, pr]] = A[[pr, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        for rr in np.nonzero(A[:, c])[0]:
            if rr != r: A[rr] = (A[rr] - A[rr, c] * A[r]) % p
        r += 1
    return r

def peval(P, al, p):
    v = 0
    for k, c in P.items():
        t = int(c) % p
        for a in k: t = t * al[a] % p
        v = (v + t) % p
    return v

def full_screen(A, ncols, al, p=2147483647 - 18):
    """is e_0 in the column span of A(al) mod p?  returns (rank A, rank [A|e0])"""
    M = np.zeros((len(A), ncols + 1), dtype=object)
    for i, R in enumerate(A):
        for c, P in R.items():
            M[i, c] = peval(P, al, p)
    M[0, ncols] = 1
    M = M.astype(np.int64) if p < 2 ** 31 else M
    return rank_mod(M[:, :ncols], p), rank_mod(M, p)
