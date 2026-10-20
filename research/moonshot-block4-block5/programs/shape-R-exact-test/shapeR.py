"""shapeR.py -- o3-shapeR (moonshot o3 round 3, 14-15 Oct 2026): solve shape R for ALL u4, u5 at once (engine game).

Setting (o3-close PROOF.md sec. 0, 3, 6). Inputs y1, y2 (any basis). Prefix u1 = D b1 y1, u2 = D b2 y2,
u3 = D(c3) with c3 in span(y, u1, u2) (shape R: extra raw terms or extra letters).  S3 = Q[B'] span(y, u1, u2, u3).
Targets: single terms T_l y_l with T_l at depth 3, one move from K, on frames F_l.  By Lemma 5+ (exactly two
T-carriers u4, u5, forced here because u1..u3 have depth <= 2): u4, u5 lie in A = S3 + Q[B']T1y1 + Q[B']T2y2 and the
2x2 target-coordinate matrix C of (u4, u5) is invertible.  Write
   M4  = { alpha : alpha.Ty in S3 + D S3 }                      (target rows of all u4 = D c4 in A)
   M5' = { alpha : alpha.Ty in S3 + D (S3 + Q[B'] T y) }        (contains the target rows of every u5, any u4 in A)
eps : Q[B'] -> Q (coefficient sum) is a ring map, so an invertible C needs eps(alpha4), eps(alpha5) independent:
   PASS(F1, F2)  <=>  eps(M4) != 0  and  dim eps(M5') = 2.     (necessary for a 2-in-7 with these targets)

THE QUOTIENT [proved in NOTE/PROOF of this folder]: Lambda = B' cap D B' D^-1 acts on the (infinite-index) closure
system by left multiplication; summing a solution over Lambda-orbits gives a solution of the finite quotient system,
and spreading a quotient solution uniformly over orbits lifts it back, so eps(M) is EXACTLY the eps-image of the
quotient.  In the engine game [B' : Lambda] = [B' : H_D] = 30 (neighbours of the input frame).  Quotient rows:
(l, F, F') with F = frame(r), F' = frame(D^-1 r) (an oriented edge of the dual polar graph); variables: S-columns
(i, N) = b_N-orbit times gen i, D-columns (i, N) = D b_N gen i, N over the 30 classes.  Entries:
   D-column (i, N): each term (l, g, q) of gen i adds q at (l, frame(D b_N g), frame(g));
   S-column (i, N): each term (l, g, q) of gen i adds q at (l, frame(g), frame(D b_N g)).
Rows with F = F_l (input l) are outputs; every other row is a constraint (= 0)."""
import sys, os, random, itertools, collections, time, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
_PL = ''
if os.path.isdir(_PL):
    sys.path.insert(0, _PL)
import numpy as np
from fractions import Fraction as Fr
import bm as BM
from fgeo import Geo

PRIMES = (1048573, 2097143)
Z4c = np.array([1, 1j, -1, -1j])
_W = np.random.default_rng(20261014).integers(1, 2**62, size=(32,), dtype=np.int64)
_W2 = np.random.default_rng(20261015).integers(1, 2**62, size=(32,), dtype=np.int64)
def vkeys(Ms):
    """(from o3-finish fulltest.py) canonical keys of B'-frames (rows up to Z4 phases and order) for a stack (n,16,16)
    on the 1/4096 grid: rows phase-normalised, hashed by two independent 64-bit linear hashes, sorted.  Equal frames
    give equal keys; different frames collide with probability ~2^-128 per pair."""
    n = Ms.shape[0]
    A = Ms.reshape(n * 16, 16)
    nz = np.abs(A) > 1e-9
    idx = np.argmax(nz, axis=1)
    v = A[np.arange(n * 16), idx]
    ang = np.angle(v)
    k = np.floor((ang + 1e-7) / (np.pi / 2)).astype(np.int64) % 4
    z = Z4c[(-k) % 4]
    A = A * z[:, None]
    re = np.round(A.real * 4096).astype(np.int64); im = np.round(A.imag * 4096).astype(np.int64)
    rows = np.concatenate([re, im], axis=1)
    with np.errstate(over='ignore'):
        h1 = (rows * _W[None, :]).sum(axis=1).reshape(n, 16)
        h2 = (rows * _W2[None, :]).sum(axis=1).reshape(n, 16)
    o = np.argsort(h1, axis=1, kind='stable')
    h1 = np.take_along_axis(h1, o, axis=1); h2 = np.take_along_axis(h2, o, axis=1)
    key = np.ascontiguousarray(np.concatenate([h1, np.sort(h2, axis=1)], axis=1))
    return list(key.view('S256').reshape(n))

def fk(M):
    return vkeys(M[None])[0]

def rankp(A, pr):
    if A.size == 0:
        return 0
    if A.shape[1] > A.shape[0]:
        A = A.T
    A = np.ascontiguousarray(np.mod(A, pr).astype(np.uint64))
    R, C = A.shape
    piv = np.zeros(max(1, min(R, C)), dtype=np.int64)
    return int(BM._lib.rref(A.ctypes.data, R, C, pr, piv.ctypes.data, 0))

def fr_mod(q, pr):
    q = Fr(q)
    return (q.numerator % pr) * pow(q.denominator % pr, pr - 2, pr) % pr

class World:
    """engine-game letters (fgeo), the 30 classes b_N, term data"""
    def __init__(self, seed=1):
        self.geo = Geo(4, seed=seed); self.G = G = self.geo.G
        self.rng = random.Random(seed)
        self.D = G.D1
        self.Di = np.conj(G.D1).T
        self.L0 = fk(G.I)
        self.Kfr = fk(G.K)
        cls = {}
        tries = 0
        while len(cls) < 30:
            b = self.geo.rand_free(); tries += 1
            k = fk(G.mul(self.D, b))
            if k not in cls:
                cls[k] = b
            assert tries < 200000
        self.classes = list(cls.values())          # b_N, N = 0..29
        self.clskey = list(cls.keys())
        self.DB = np.stack([G.mul(self.D, b) for b in self.classes])     # D b_N, exact
        self.cache = {}
        self.dep = {}

    def st(self, M):
        return self.geo.d0(M), self.geo.dK(M)

    def nbrs(self, M):
        """frames of D b_N M for the 30 classes (memoised by matrix key)"""
        G = self.G
        k = G.key(M)
        v = self.cache.get(k)
        if v is None:
            P = np.matmul(self.DB, M[None])
            P = np.round(P.real * 4096) / 4096 + 1j * (np.round(P.imag * 4096) / 4096)
            v = (fk(M), vkeys(P))
            self.cache[k] = v
            if v[0] not in self.dep:
                self.dep[v[0]] = float(np.log2((np.abs(M) > 1e-9).sum(axis=1).max()))
            rs = (np.abs(P) > 1e-9).sum(axis=2).max(axis=1)
            for kk, d in zip(v[1], rs):
                if kk not in self.dep:
                    self.dep[kk] = float(np.log2(d))
        return v

    def rand_letter(self):
        return self.geo.rand_free()

# generators: list of terms (l, M, q) with M a 16x16 matrix (exact grid), q a Fraction
def gen_inp(l):
    return [(l, None, Fr(1))]

def lmul(G, b, gen, I):
    return [(l, G.mul(b, I if M is None else M), q) for (l, M, q) in gen]

def collect(G, terms, I):
    """merge equal matrices"""
    out = {}
    mats = {}
    for (l, M, q) in terms:
        M = I if M is None else M
        k = (l, G.key(M))
        out[k] = out.get(k, 0) + Fr(q); mats[k] = M
    return [(k[0], mats[k], v) for k, v in out.items() if v != 0]

def move(W, b, content):
    """u = D b (content); content = list of terms"""
    G = W.G
    return collect(G, [(l, G.mul(W.D, b, (G.I if M is None else M)), q) for (l, M, q) in content], G.I)

def comb(*pairs):
    out = []
    for (q, gen) in pairs:
        out += [(l, M, Fr(q) * c) for (l, M, c) in gen]
    return out

def eps_system(W, Sgens, Dgens, Fl, lev=0):
    """quotient system; returns (constraint rows, output rows per input l) as sparse dicts var -> coef.
    lev > 0: RELAXED level -- constraint rows on frames of depth < lev are dropped (junk below depth lev ignored)."""
    G = W.G
    rows = collections.defaultdict(lambda: collections.defaultdict(Fr))
    var = {}
    def v(key):
        if key not in var:
            var[key] = len(var)
        return var[key]
    for i, gen in enumerate(Sgens):
        for (l, M, q) in gen:
            fr, nb = W.nbrs(G.I if M is None else M)
            for N in range(30):
                rows[(l, fr, nb[N])][v(('S', i, N))] += q
    for i, gen in enumerate(Dgens):
        for (l, M, q) in gen:
            fr, nb = W.nbrs(G.I if M is None else M)
            for N in range(30):
                rows[(l, nb[N], fr)][v(('D', i, N))] += q
    cons = []; outs = {l: collections.defaultdict(Fr) for l in Fl}
    for (l, F, F2), r in rows.items():
        if F == Fl.get(l):
            for c, q in r.items():
                outs[l][c] += q
        elif lev > 0 and W.dep[F] < lev - 1e-9:
            continue
        else:
            rr = {c: q for c, q in r.items() if q != 0}
            if rr:
                cons.append(rr)
    return cons, outs, len(var)

def eps_dim(cons, outs, nv, primes=PRIMES):
    """dim eps(M) = rank[C; e_1..e_W] - rank C (max over primes), then the same for each single output e_l"""
    res = []
    L = sorted(outs)
    for pr in primes:
        A = np.zeros((len(cons) + len(L), nv), dtype=np.float64)
        for i, r in enumerate(cons):
            for c, q in r.items():
                A[i, c] = fr_mod(q, pr)
        for j, l in enumerate(L):
            for c, q in outs[l].items():
                A[len(cons) + j, c] = fr_mod(q, pr)
        rC = rankp(A[:len(cons)], pr) if cons else 0
        rall = rankp(A, pr)
        singles = []
        for j in range(len(L)):
            B = np.vstack([A[:len(cons)], A[len(cons) + j:len(cons) + j + 1]])
            singles.append(rankp(B, pr) - rC)
        res.append(tuple([rall - rC] + singles))
    return max(res)

def target_frames(W, gens, l):
    """frames F at depth 3, one move from K (dist 1), of the form D b_N g for a depth-2 term g of gens on input l"""
    G = W.G
    out = {}
    for gen in gens:
        for (ll, M, q) in gen:
            if ll != l or M is None:
                continue
            if W.st(M)[0] != 2:
                continue
            for N, b in enumerate(W.classes):
                T = G.mul(W.D, b, M)
                if W.st(T) == (3, 1):
                    k = fk(T)
                    out.setdefault(k, T)
    return out

def test_prefix(W, gens3, verbose=False, lev=0):
    """gens3 = [y1, y2, u1, u2, u3] (S3 generators). Returns per target pair (F1, F2): (dimeps M4, dimeps M5')"""
    T1s = target_frames(W, gens3, 0); T2s = target_frames(W, gens3, 1)
    res = []
    for (k1, T1), (k2, T2) in itertools.product(T1s.items(), T2s.items()):
        Fl = {0: k1, 1: k2}
        c4, o4, n4 = eps_system(W, gens3, gens3, Fl, lev)
        d4 = eps_dim(c4, o4, n4)
        Ty = [[(0, T1, Fr(1))], [(1, T2, Fr(1))]]
        c5, o5, n5 = eps_system(W, gens3, gens3 + Ty, Fl, lev)
        d5 = eps_dim(c5, o5, n5)
        res.append({"d4": d4, "d5": d5, "pass": d4[0] >= 1 and d5[0] == 2})
    return res, len(T1s), len(T2s)
