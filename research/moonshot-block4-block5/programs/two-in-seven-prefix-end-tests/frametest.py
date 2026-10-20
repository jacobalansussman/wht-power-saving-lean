"""frametest.py -- o3-finish2 (14 Oct 2026): the FRAME-LEVEL (augmentation) relaxation of the 2-in-7 end test.
ONE finite linear system per prefix and per one-dimensional character chi of B' bounds, for EVERY end move E and every
invertible mixing Q in GL_2(Q[B']) at once, the augmentation of the end rows (NOTE.md of o3-finish2, Lemma F).

Setting (widened game, as bm.py / o3lib.py): Gamma = invertible 16x16 matrices on the dyadic grid, B' = Z4-monomials,
frames = left cosets B'g.  Fix one canonical representative f_F per frame (rows phase-normalised, sorted).  For a
character chi of B' (chi_{a,k}(b) = sgn(perm b)^a * (prod of b's entries)^k, a in {0,1}, k in Z4: all 8 linear
characters of Z4 wr S16) define the left-Q[B']-semilinear map
        pi_chi : Q[Gamma] y1 + Q[Gamma] y2 -> Q(i)[frames x {1,2}],   b f_F y_l |-> chi(b) [F, l].
Then pi_chi(q z) = chi(q) pi_chi(z) for q in Q[B'].  For a word with generators y1, y2, u1..u5 and an end
        a1 K y1 + a2 K y2 = s + E s',   s, s' in S = sum Q[B'] gens,  E = c D1 c' any paid element,
applying pi_chi gives chi(a1) kappa_1 + chi(a2) kappa_2 (kappa_l = chi(K f_K^-1)[K-frame, l]) in the Q(i)-span V_chi of
  * span columns  pi_chi(u_t)                         (t = 1..5; raw inputs only touch the free frame),
  * D1 columns    pi_chi(D1 gamma u_t), gamma in B'  (these depend on gamma only through the coset H_D gamma, up to a
                  scalar; for a row frame F and a term g of u_t with h = f_F g^-1 one move, D1 gamma = b h, so the column
                  is pi_chi(h u_t) up to the scalar chi(b)),
  * rows on the free frame and on frames one move from it are covered by single-entry columns (raw inputs, D1 gamma y_l)
    and are dropped exactly.
So E_chi := { e in Q(i)^2 : e1 kappa_1 + e2 kappa_2 in V_chi } contains chi(M_E(B')) for EVERY E.  A 2-in-7 with end
K x = Q (E6 c6, E7 c7)^T + s, Q invertible over Q[B'], has the rows of Q^-1 in M_E6, M_E7, and chi(Q^-1) invertible
(chi is a ring map), so it needs dim E_chi = 2 for every chi.  dim E_chi <= 1 for ONE chi excludes the prefix for every
end move, every basis and every mixing.

The closure is explored breadth first from the two K-frame rows; at a cap the relaxation (processed rows, every column
touching them) is used (sound: fewer equations).  Linear algebra mod p (two primes, i -> sqrt(-1) mod p) for screening;
exact certificates (left null vectors over Q(i), verified exactly) in frameexact.py.
usage: python3 frametest.py PREFIXES.jsonl OUT.jsonl [first] [step] [caprows] [only RUN:IT,RUN:IT,...]"""
import sys, os, json, time, collections
PYLIB = ''
if os.path.isdir(PYLIB):
    sys.path.insert(0, PYLIB)
HERE = os.path.dirname(os.path.abspath(__file__))
V3 = os.environ.get('O3_V3', os.path.join(HERE, 'v3'))
if not os.path.isdir(V3):
    raise SystemExit('v3/ not found next to frametest.py')
sys.path.insert(0, V3)
import numpy as np
from fractions import Fraction
from o3lib import W4, build
import bm as BM

SC = 4096
PRIMES = (1048573, 2097133)


def sqrtm1(p):
    for g in range(2, 5000):
        if pow(g, (p - 1) // 2, p) == p - 1:
            return pow(g, (p - 1) // 4, p)
    raise ValueError


# ------------------------------------------------------------------ integer matrices
def toint(M):
    re = np.round(M.real * SC).astype(np.int64); im = np.round(M.imag * SC).astype(np.int64)
    assert np.abs(re - M.real * SC).max() < 1e-6 and np.abs(im - M.imag * SC).max() < 1e-6
    return re, im


def frame_data(M):
    """canonical frame of B'M: (key, f, parity, rsum) with M = b f, b monomial, chi_{a,k}(b) = (-1)^(a parity) i^(-k rsum)"""
    re, im = toint(M)
    n = re.shape[0]
    r = np.zeros(n, dtype=np.int64)
    R = re.copy(); I = im.copy()
    for x in range(n):
        nz = np.nonzero((R[x] != 0) | (I[x] != 0))[0][0]
        a, b = R[x, nz], I[x, nz]
        k = 0
        while not (a > 0 and b >= 0):
            a, b = -b, a          # multiply by i
            k += 1
        # apply i^k to the whole row
        for _ in range(k):
            R[x], I[x] = -I[x].copy(), R[x].copy()
        r[x] = k
    rows = np.concatenate([R, I], axis=1)
    order = sorted(range(n), key=lambda x: tuple(rows[x]))          # f_j = normalised row order[j]
    key = rows[order].tobytes()
    # b[x, pi^-1(x)] = i^(-r_x): permutation sigma = pi where f_j = n_{pi(j)} -> pi(j) = order[j]
    perm = order
    seen = [False] * n; par = 0
    for s in range(n):
        if not seen[s]:
            L = 0; j = s
            while not seen[j]:
                seen[j] = True; j = perm[j]; L += 1
            par ^= (L - 1) & 1
    f = (R[order] + 1j * I[order]) / SC
    return key, f, par, int(r.sum() % 4)


def is_mono_int(re, im):
    nz = (re != 0) | (im != 0)
    if not ((nz.sum(axis=1) == 1).all() and (nz.sum(axis=0) == 1).all()):
        return False
    v = np.abs(re[nz]) + np.abs(im[nz])
    return bool((v == SC).all() and ((re[nz] == 0) | (im[nz] == 0)).all())


def one_move_int(re, im):
    """h in B' D1 B' ?  (rows support 2, entries in A Z4, rows pair up by support, supports partition the columns,
    ratio of the partner row = - ratio)"""
    nz = (re != 0) | (im != 0)
    if not (nz.sum(axis=1) == 2).all():
        return False
    if not ((np.abs(re[nz]) == SC // 2).all() and (np.abs(im[nz]) == SC // 2).all()):
        return False
    sup = {}
    for x in range(re.shape[0]):
        c = tuple(np.nonzero(nz[x])[0])
        sup.setdefault(c, []).append(x)
    cols = set()
    for c, xs in sup.items():
        if len(xs) != 2:
            return False
        cols.update(c)
    if len(cols) != re.shape[1]:
        return False
    def m(a, b):        # 2048 (1+i) i^m
        return {(1, 1): 0, (-1, 1): 1, (-1, -1): 2, (1, -1): 3}[(int(np.sign(a)), int(np.sign(b)))]
    for (c1, c2), (x, x2) in sup.items():
        rho = (m(re[x, c2], im[x, c2]) - m(re[x, c1], im[x, c1])) % 4
        rho2 = (m(re[x2, c2], im[x2, c2]) - m(re[x2, c1], im[x2, c1])) % 4
        if (rho2 - rho) % 4 != 2:
            return False
    return True


class FrameSys:
    def __init__(self, rec, caprows=200000):
        self.w = W4(1); G = self.G = self.w.G
        self.gens = build(self.w, rec["prog"])
        self.caprows = caprows
        self.fd = {}          # bytes of matrix -> frame data (+ freeish flag)
        self.fmat = {}        # frame key -> f
        self.terms = []       # per gen t: list of (l, M, coef)
        for t, g in enumerate(self.gens):
            self.terms.append([(l, G.MATS[gid], Fraction(c)) for (l, gid), c in g.items()])
        self.index = collections.defaultdict(list)    # (key, l) -> [(t, j)]
        for t in range(2, len(self.gens)):
            for j, (l, M, c) in enumerate(self.terms[t]):
                d = self.data(M)
                self.index[(d[0], l)].append((t, j))
        self.Kkey = self.data(G.K)[0]

    def data(self, M):
        kb = toint(M)
        bkey = kb[0].tobytes() + kb[1].tobytes()
        d = self.fd.get(bkey)
        if d is None:
            re, im = kb
            freeish = is_mono_int(re, im) or one_move_int(re, im)
            key, f, par, rs = frame_data(M)
            self.fmat.setdefault(key, f)
            d = (key, par, rs, freeish)
            self.fd[bkey] = d
        return d

    def column(self, t, h):
        """entries of pi(h u_t) (h None: pi(u_t)): {(key, l): [(coef, par, rsum), ...]} without freeish rows"""
        G = self.G
        out = collections.defaultdict(list)
        for (l, M, c) in self.terms[t]:
            Mh = M if h is None else G.mul(h, M)
            key, par, rs, freeish = self.data(Mh)
            if freeish:
                continue
            out[(key, l)].append((c, par, rs))
        return out

    def run(self):
        G = self.G
        rows = {}            # (key, l) -> index
        cols = {}            # column id -> entries
        q = collections.deque(); processed = []
        def addrow(rk):
            if rk not in rows:
                rows[rk] = len(rows); q.append(rk)
        for l in range(2):
            addrow((self.Kkey, l))
        capped = False
        while q:
            if len(processed) >= self.caprows:
                capped = True; break
            rk = q.popleft()
            processed.append(rk)
            key, l = rk
            f = self.fmat[key]
            # span columns
            for (t, j) in self.index.get(rk, ()):
                cid = ('S', t)
                if cid not in cols:
                    cols[cid] = self.column(t, None)
                    for r2 in cols[cid]:
                        addrow(r2)
            # D1 columns through every term of every move on input l
            for t in range(2, len(self.gens)):
                for j, (lj, M, c) in enumerate(self.terms[t]):
                    if lj != l:
                        continue
                    h = G.mul(f, np.conj(M).T)
                    re, im = toint(h)
                    if not one_move_int(re, im):
                        continue
                    hk = self.data(h)[0]
                    cid = ('D', t, hk)
                    if cid in cols:
                        continue
                    cols[cid] = self.column(t, h)
                    for r2 in cols[cid]:
                        addrow(r2)
        self.rows, self.cols, self.processed, self.capped = rows, cols, processed, capped
        return self

    def matrix_modp(self, a, k, p, ip):
        """entries mod p for chi_{a,k}; rows = processed rows only (relaxation if capped); 2 seed columns last"""
        keep = {rk: i for i, rk in enumerate(self.processed)}
        ent = collections.defaultdict(int)
        cl = list(self.cols.keys())
        ipow = [1, ip, ip * ip % p, ip * ip * ip % p]
        for ci, cid in enumerate(cl):
            for rk, lst in self.cols[cid].items():
                ri = keep.get(rk)
                if ri is None:
                    continue
                v = 0
                for (c, par, rs) in lst:
                    s = (c.numerator % p) * pow(c.denominator % p, p - 2, p) % p
                    if a and par:
                        s = (-s) % p
                    s = s * ipow[(-k * rs) % 4] % p
                    v = (v + s) % p
                if v:
                    ent[(ri, ci)] = v
        nc = len(cl)
        for l in range(2):
            ent[(keep[(self.Kkey, l)], nc + l)] = 1
        return len(keep), nc, ent


def dimE_modp(nr, nc, ent, p):
    """dim E = 2 - (rank[A S] - rank A), with exact singleton peeling first (as fulltest.solve_dimM)"""
    rows = collections.defaultdict(dict); colsd = collections.defaultdict(set)
    for (r, c), v in ent.items():
        rows[r][c] = v; colsd[c].add(r)
    T = {nc, nc + 1}
    free_T = sum(1 for t in T if not colsd.get(t))
    T = {t for t in T if colsd.get(t)}
    changed = True
    while changed:
        changed = False
        for r in list(rows.keys()):
            rr = rows.get(r)
            if rr is None:
                continue
            if len(rr) == 0:
                del rows[r]; changed = True; continue
            if len(rr) == 1:
                (c,) = rr.keys()
                for r2 in list(colsd[c]):
                    rows[r2].pop(c, None)
                colsd[c] = set(); T.discard(c)
                del rows[r]; changed = True
        for c in list(colsd.keys()):
            if c in T:
                continue
            if len(colsd[c]) == 1:
                (r,) = colsd[c]
                if r in rows:
                    for c2 in list(rows[r].keys()):
                        colsd[c2].discard(r)
                    del rows[r]
                colsd[c] = set(); changed = True
    free_T += sum(1 for t in T if not colsd.get(t))
    T = {t for t in T if colsd.get(t)}
    rl = list(rows.keys())
    if not T:
        return free_T, (0, 0)
    cl = sorted({c for r in rl for c in rows[r]})
    Xc = [c for c in cl if c not in T]; Tc = [c for c in cl if c in T]
    if len(rl) * len(cl) > 4e7:
        return None, (len(rl), len(cl))
    A = np.zeros((len(rl), len(Tc) + len(Xc)), dtype=np.float64)
    ci = {c: k for k, c in enumerate(Tc + Xc)}
    for i, r in enumerate(rl):
        for c, v in rows[r].items():
            A[i, ci[c]] = v
    rA = rankp(A, p) if A.size else 0
    rX = rankp(A[:, len(Tc):], p) if len(Xc) else 0
    return free_T + len(Tc) - rA + rX, (len(rl), len(cl))


def rankp(A, pr):
    if A.shape[1] > A.shape[0]:
        A = A.T
    A = np.ascontiguousarray(np.mod(A, pr).astype(np.uint64))
    R, C = A.shape
    piv = np.zeros(max(1, min(R, C)), dtype=np.int64)
    return int(BM._lib.rref(A.ctypes.data, R, C, pr, piv.ctypes.data, 0))


CHARS = [(a, k) for a in (0, 1) for k in range(4)]


def test(rec, caprows=200000):
    t0 = time.time()
    fs = FrameSys(rec, caprows).run()
    out = {"run": rec.get("run", "").replace(".jsonl", ""), "it": rec.get("it"), "rows": len(fs.rows),
           "proc": len(fs.processed), "cols": len(fs.cols), "capped": fs.capped, "dims": {}}
    for (a, k) in CHARS:
        ds = []
        for p in PRIMES:
            ip = sqrtm1(p)
            nr, nc, ent = fs.matrix_modp(a, k, p, ip)
            d, core = dimE_modp(nr, nc, ent, p)
            ds.append(d)
        out["dims"]["%d%d" % (a, k)] = ds
        out.setdefault("core", core)
    mins = [max(v) if None not in v else None for v in out["dims"].values()]
    out["min_dim"] = min([m for m in mins if m is not None], default=None)
    out["sec"] = round(time.time() - t0, 1)
    return out


if __name__ == "__main__":
    assert all(p % 4 == 1 for p in PRIMES), PRIMES
    inp, outp = sys.argv[1], sys.argv[2]
    first = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    step = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    caprows = int(sys.argv[5]) if len(sys.argv) > 5 else 200000
    only = None
    if len(sys.argv) > 6:
        only = {(x.split(':')[0], int(x.split(':')[1])) for x in sys.argv[6].split(',')}
    recs = [json.loads(l) for l in open(inp)]
    fo = open(outp, "a")
    for kk in range(first, len(recs), step):
        rec = recs[kk]
        key = (rec.get("run", "").replace(".jsonl", ""), rec.get("it"))
        if only is not None and key not in only:
            continue
        r = test(rec, caprows)
        r["k"] = kk
        fo.write(json.dumps(r) + "\n"); fo.flush()
        print(kk, r["run"], r["it"], "rows", r["rows"], "proc", r["proc"], "cols", r["cols"], "capped", r["capped"],
              "dims", r["dims"], "min", r["min_dim"], r["sec"], "s", flush=True)
