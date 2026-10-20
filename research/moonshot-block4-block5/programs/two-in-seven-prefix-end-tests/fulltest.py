"""fulltest.py -- o3-finish (14 Oct 2026): the COMPLETE any-basis private-end test for a saved prefix.

For a prefix (program record from o3f-replay / o3f-cov2), a paid element E (or E = None: no end move) and the basis
group H = P = S4 x| Pauli4 (bit permutations, X_a, Z_b and the scalar letter iI: 24,576 letters, all normalising K;
it contains R = the 384 address maps of o3-build/anybasis.py), compute
    M_E(H) = { (a, b) in Q[H]^2 : K (a x1 + b x2) in span + E span }.
Two private ends in SOME basis y = F x, F in GL_2(Q[H]), with end moves E6, E7 need M_E6 != 0 and M_E7 != 0. This
covers PRIV (F = I), the +-1 halving end and the order-3 halving end for every pair (E6, E7) (o3-build Cor. 1b), so
"M_E(H) = 0 for every class of E" excludes all of them for that prefix at once.

Exactness and completeness (no cap):
  * E runs over ALL paid elements of the prefix (search.E_pool without the limit 3), one per class modulo
    E ~ lam E b (lam in H, b free): M_{lam E} = mu M_E with mu in H, so classes are equivalent.  E = None is added.
  * Seeds (l, s K), s in H.  A seed row no column touches forces its coefficient to 0 (dropped, exact).  The others
    split into orbits of Lambda_E = { lam in H : E^-1 lam E free } acting on the left; Lambda_E maps the closure
    system to itself, so the system is a direct sum of connected components that Lambda_E permutes, and
    M_E = 0 iff the component of ONE seed per orbit has M = 0.
  * Each component is explored breadth first (rows nearest the seeds first, seed rows always first).  If it closes
    below the cap, the system is exact.  If not, the RELAXATION (processed rows only, every column that touches them)
    is solved: a relaxation with M = 0 proves M = 0 (fewer equations can only enlarge M).  Only a relaxation with
    M != 0 at the largest cap is reported as undecided.
  * Linear algebra mod p = 1048573 (bm.P): singleton peeling (exact), then dense rank of the core.  dim M =
    #T - rank[A_T A_X] + rank A_X.  A nonzero answer is re-checked with a second prime before it is reported.
Free rows (monomial) and E-free rows (E^-1 g monomial) are dropped exactly: each is covered by a single-term column on
an input that touches no other row.
usage: python3 -I fulltest.py PREFIXES.jsonl OUT.jsonl [first] [step] [caprows]"""
import sys, os, json, time, itertools, collections, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'v3'))
from o3lib import *
import bm as BM
import search as S

P_ = BM.P
Z4c = np.array([1, 1j, -1, -1j])

# ------------------------------------------------------------------ monomials as (perm, phase) arrays
def mono_pq(G, M):
    s, q = G.mono_of(M)
    return np.array(s, dtype=np.int64), np.array(q, dtype=np.int64) % 4

def pq_mat(p, q):
    M = np.zeros((16, 16), dtype=complex)
    M[p, np.arange(16)] = Z4c[q % 4]
    return M

def H_group():
    """P = S4 x| Pauli4 with phases: e_y -> i^c (-1)^{b.sig(y)} e_{sig(y)^a}; returns arrays p (n,16), q (n,16), and
    flags isR (c = 0, b = 0: the 384 address maps of anybasis.py)"""
    ps, qs, isR = [], [], []
    for perm in itertools.permutations(range(4)):
        sig = perm_bits(list(perm))
        for a in range(16):
            for b in range(16):
                for c in range(4):
                    ps.append([sig[y] ^ a for y in range(16)])
                    qs.append([(c + 2 * dot(b, sig[y])) % 4 for y in range(16)])
                    isR.append(b == 0 and c == 0)
    return np.array(ps), np.array(qs), np.array(isR)

def in_H(p, q):
    """monomial (p, q) in P ?"""
    a = int(p[0]); sig = p ^ a
    imgs = [int(sig[1 << j]) for j in range(4)]
    if sorted(imgs) != [1, 2, 4, 8]:
        return False
    for y in range(16):
        v = 0
        for j in range(4):
            if (y >> j) & 1: v ^= imgs[j]
        if v != sig[y]:
            return False
    d = (q - q[0]) % 4
    if np.any(d % 2):
        return False
    bits = (d // 2) % 2          # must be b . sig(y)
    bvec = 0
    inv = {int(sig[y]): y for y in range(16)}
    for j in range(4):
        if bits[inv[1 << j]]: bvec |= 1 << j
    return all(int(bits[y]) == dot(bvec, int(sig[y])) for y in range(16))

# ------------------------------------------------------------------ vectorised canonical frame keys (own convention)
_W = np.random.default_rng(20261014).integers(1, 2**62, size=(32,), dtype=np.int64)
_W2 = np.random.default_rng(20261015).integers(1, 2**62, size=(32,), dtype=np.int64)
def vkeys(Ms):
    """canonical keys of B'-frames (left monomial orbits: rows up to Z4 phases and order) for a stack (n,16,16) of
    matrices on the 1/4096 grid: each row is phase-normalised (first nonzero entry's angle in [0, pi/2)) and hashed
    (two independent 64-bit linear hashes, wrapping arithmetic), the 16 row hashes are sorted.  Equal frames give equal
    keys deterministically; different frames collide with probability ~2^-128 per pair."""
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

def is_mono_stack(Ms):
    nz = np.abs(Ms) > 1e-9
    ok = (nz.sum(axis=2) == 1).all(axis=1) & (nz.sum(axis=1) == 1).all(axis=1)
    vals = np.where(nz, Ms, 0).sum(axis=2)
    ok &= (np.abs(np.abs(vals) - 1) < 1e-9).all(axis=1)
    return ok

# ------------------------------------------------------------------ sparse rank with peeling
def rankp(A, pr=BM.P):
    """rank mod pr (pr < 2^21: bm's C rref uses delayed reduction)"""
    if A.shape[1] > A.shape[0]:
        A = A.T
    A = np.ascontiguousarray(np.mod(A, pr).astype(np.uint64))
    R, C = A.shape
    piv = np.zeros(max(1, min(R, C)), dtype=np.int64)
    return int(BM._lib.rref(A.ctypes.data, R, C, pr, piv.ctypes.data, 0))

PRIMES = tuple(int(x) for x in os.environ.get('O3F_PRIMES', str(BM.P)).split(','))

def solve_dimM(nrows, ent, tset, primes=None):
    """ent: list of (row, col, Fraction); tset: set of T columns. Returns dim M = #T - rank(A) + rank(A_X) (mod p),
    after exact singleton peeling.  'free T' (a T column in no row) counts as +1 directly."""
    res = []
    if primes is None:
        primes = PRIMES
    for pr in primes:
        rows = collections.defaultdict(dict); cols = collections.defaultdict(set)
        for (r, c, co) in ent:
            co = Fraction(co)
            v = (co.numerator % pr) * pow(co.denominator % pr, pr - 2, pr) % pr
            nv = (rows[r].get(c, 0) + v) % pr
            if nv: rows[r][c] = nv; cols[c].add(r)
            else:
                rows[r].pop(c, None); cols[c].discard(r)
        T = set(tset)
        free_T = sum(1 for t in T if not cols.get(t))
        T = {t for t in T if cols.get(t)}
        # peeling
        changed = True
        while changed:
            changed = False
            for r in list(rows.keys()):
                rr = rows.get(r)
                if rr is None: continue
                if len(rr) == 0:
                    del rows[r]; changed = True; continue
                if len(rr) == 1:
                    (c,) = rr.keys()
                    # forces column c = 0: remove column c everywhere
                    for r2 in list(cols[c]):
                        rows[r2].pop(c, None)
                    cols[c] = set(); T.discard(c)
                    del rows[r]; changed = True
            for c in list(cols.keys()):
                if c in T: continue
                if len(cols[c]) == 1:
                    (r,) = cols[c]
                    if r in rows:
                        for c2 in list(rows[r].keys()):
                            cols[c2].discard(r)
                        del rows[r]
                    cols[c] = set(); changed = True
        free_T += sum(1 for t in T if not cols.get(t))
        T = {t for t in T if cols.get(t)}
        rl = list(rows.keys())
        if not T:
            res.append(free_T); continue
        cl = sorted({c for r in rl for c in rows[r]})
        Xc = [c for c in cl if c not in T]; Tc = [c for c in cl if c in T]
        if len(rl) * len(cl) > MAXCORE:
            res.append(None); continue
        A = np.zeros((len(rl), len(Tc) + len(Xc)), dtype=np.float64)
        ci = {c: k for k, c in enumerate(Tc + Xc)}
        for i, r in enumerate(rl):
            for c, v in rows[r].items():
                A[i, ci[c]] = v
        rA = rankp(A, pr) if A.size else 0
        rX = rankp(A[:, len(Tc):], pr) if len(Xc) else 0
        res.append(free_T + len(Tc) - rA + rX)
    if any(x is None for x in res):
        return None
    return max(res)          # several primes: zero only if zero for every prime

# ------------------------------------------------------------------ BFS closure for one component
class Ctx:
    def __init__(self, w, gens):
        self.w = w; self.G = w.G; self.gens = gens
        G = self.G
        self.index = {}
        for i, g in enumerate(gens):
            for (ll, gid) in g:
                self.index.setdefault((G.frm(gid), ll), []).append((i, gid))
        self.Kfrm = G.frm(G.KID)
        self.Ki = np.conj(G.K).T
        self._seedmemo = {}

    def seed_of(self, gid):
        """row g is a seed (g = s K, s in P)?"""
        v = self._seedmemo.get(gid)
        if v is None:
            G = self.G
            v = False
            if G.frm(gid) == self.Kfrm:
                m = G.mul(G.MATS[gid], self.Ki)
                if G.is_monomial(m):
                    p, q = mono_pq(G, m)
                    v = in_H(p, q) and (not ONLY_R or not np.any(q % 4))
            self._seedmemo[gid] = v
        return v

CHECKS = (300, 1500, 6000, 20000)
MATCAP = int(os.environ.get('O3F_MATCAP', 250000))   # interned matrices (~12 KB each): memory guard
MAXCORE = float(os.environ.get('O3F_MAXCORE', 4e6))   # dense core after peeling: larger -> undecided at this checkpoint

def component(ctx, E, start, caprows, capcols):
    """BFS closure of the homogeneous seed system from the seed row `start` = (l, gid) (an orbit representative).
    Verdict 'zero' iff the representative's own coordinate t_start is forced to 0, in the exact component or in a
    relaxation (processed rows + every column touching them) at a checkpoint: by the Lambda_E-symmetry every seed is a
    translate of its orbit's representative, so all representatives forced => M_E = 0 (NOTE.md, section 2).
    'nonzero': exact component with t_start free (M_E != 0).  'undecided': relaxation at the cap with t_start free."""
    G = ctx.G; gens = ctx.gens; index = ctx.index
    Einv = None if E is None else np.conj(E).T
    rows = {}; cols = {}; ent = []; tset = set()
    q = collections.deque(); processed = set()
    def freeish(g3):
        if G.is_free(g3):
            return True
        return Einv is not None and G.is_monomial(G.mul(Einv, G.MATS[g3]))
    def addrow(t):
        r = rows.get(t)
        if r is None:
            r = len(rows); rows[t] = r
            if ctx.seed_of(t[1]):
                c = ('T', t); cols[c] = len(cols); tset.add(cols[c]); ent.append((r, cols[c], -1))
                q.appendleft(t)
            else:
                q.append(t)
        return r
    addrow(start)
    tc = cols[('T', start)]
    complete = True
    checks = [c for c in CHECKS if c < caprows] + [caprows]
    ci = 0; last = None
    while q:
        if ci < len(checks) and len(processed) >= checks[ci]:
            ci += 1
            keep = {rows[t] for t in processed}
            use = [(r, c, co) for (r, c, co) in ent if r in keep]
            d = solve_dimM(len(rows), use, {tc})
            last = d
            if d == 0:
                return {"status": "relaxed", "verdict": "zero", "rows": len(rows), "proc": len(processed),
                        "cols": len(cols), "nT": len(tset), "dimM": 0}
            if ci == len(checks):
                complete = False; break
        if len(cols) >= capcols or len(G.MATS) > MATCAP:
            complete = False; break
        t = q.popleft()
        if t in processed: continue
        processed.add(t)
        (lt, gt) = t
        Mg = G.MATS[gt]
        cands = []
        for (i, gterm) in index.get((G.frm(gt), lt), ()):
            cands.append(('S', i, G.mul(Mg, np.conj(G.MATS[gterm]).T), None))
        if E is not None:
            hm = G.mul(Einv, Mg)
            for (i, gterm) in index.get((G.frame_key(hm), lt), ()):
                cands.append(('E', i, G.mul(hm, np.conj(G.MATS[gterm]).T), E))
        for (kind, i, bmat, pre) in cands:
            if not G.is_monomial(bmat):
                continue
            b = G.intern(bmat); ck = (kind, i, b)
            if ck in cols:
                continue
            c = len(cols); cols[ck] = c
            P2 = bmat if pre is None else G.mul(pre, bmat)
            for (l2, g2), co in gens[i].items():
                g3 = G.intern(G.mul(P2, G.MATS[g2]))
                if freeish(g3):
                    continue
                ent.append((addrow((l2, g3)), c, co))
    if complete:
        d = solve_dimM(len(rows), ent, {tc})
        dall = solve_dimM(len(rows), ent, set(tset))
        verdict = "zero" if d == 0 else ("nonzero" if d is not None else "undecided")
        return {"status": "exact", "verdict": verdict, "rows": len(rows), "proc": len(processed),
                "cols": len(cols), "nT": len(tset), "dimM": dall, "dim_start": d}
    if last is None:
        keep = {rows[t] for t in processed}
        last = solve_dimM(len(rows), [(r, c, co) for (r, c, co) in ent if r in keep], {tc})
    return {"status": "relaxed", "verdict": "zero" if last == 0 else "undecided", "rows": len(rows),
            "proc": len(processed), "cols": len(cols), "nT": len(tset), "dimM": last}

# ------------------------------------------------------------------ per-prefix driver
def e_classes(w, gens, Hp, Hq, HisR):
    """all paid elements (E_pool unlimited) with their R-class and P-class keys"""
    G = w.G
    Es = S.E_pool(w, gens, 10 ** 9)
    Hm = np.zeros((len(Hp), 16, 16), dtype=complex)
    Hm[np.arange(len(Hp))[:, None], Hp, np.arange(16)[None, :]] = Z4c[Hq]
    Hinv = np.conj(np.transpose(Hm, (0, 2, 1)))
    out = []
    for E in Es:
        Ei = np.conj(E).T
        X = np.matmul(Ei[None], Hinv)        # E^-1 lam^-1 for all lam
        ks = vkeys(X)
        kR = min(k for k, f in zip(ks, HisR) if f)
        kP = min(ks)
        out.append((E, kR, kP))
    return out

def seed_orbits(w, ctx, E, l, Hp, Hq, Hm, Hkeyidx):
    """alive seeds s K (s in P) of array l and one representative per Lambda_E orbit"""
    G = w.G
    n = len(Hp)
    SK = np.matmul(Hm, G.K[None])
    alive = np.zeros(n, dtype=bool)
    if (ctx.Kfrm, l) in ctx.index:
        alive[:] = True
    if E is not None:
        Ei = np.conj(E).T
        ks = vkeys(np.matmul(Ei[None], SK))
        tk = set()
        for i, g in enumerate(ctx.gens):
            for (ll, gid) in g:
                if ll == l:
                    tk.add(vkeys(G.MATS[gid][None])[0])
        alive |= np.array([k in tk for k in ks])
    if E is None:
        lam = np.arange(n)
    else:
        Ei = np.conj(E).T
        conj = np.matmul(np.matmul(Ei[None], Hm), E[None])
        lam = np.nonzero(is_mono_stack(conj))[0]
    # orbits: s -> lam s
    seen = np.zeros(n, dtype=bool); reps = []
    for s in np.nonzero(alive)[0]:
        if seen[s]: continue
        reps.append(int(s))
        # lam s as (p, q): p = p_lam[p_s], q = q_s + q_lam[p_s]
        ps = Hp[s]; qs = Hq[s]
        P2 = Hp[lam][:, ps]; Q2 = (qs[None, :] + Hq[lam][:, ps]) % 4
        for k in range(len(lam)):
            j = Hkeyidx.get((P2[k].tobytes(), Q2[k].tobytes()))
            if j is not None: seen[j] = True
    return int(alive.sum()), len(lam), reps, SK

def test_prefix(rec, caprows=60000, capcols=60000, Hdata=None):
    w = W4(1); G = w.G
    gens = build(w, rec["prog"])
    ctx = Ctx(w, gens)
    Hp, Hq, HisR, Hm, Hkeyidx = Hdata
    t0 = time.time()
    ecl = e_classes(w, gens, Hp, Hq, HisR)
    Rcls = sorted({kR for (_, kR, _) in ecl}); Pcls = {}
    for (E, kR, kP) in ecl:
        Pcls.setdefault(kP, E)
    tests = [("None", None)] + [("E%d" % k, E) for k, E in enumerate(Pcls.values())]
    out = {"run": rec.get("run"), "it": rec.get("it"), "nE": len(ecl), "nRclass": len(Rcls), "nPclass": len(Pcls),
           "tests": []}
    allzero = True; undecided = False
    for (nm, E) in tests:
        trec = {"E": nm, "orbits": []}
        for l in range(2):
            na, nlam, reps, SK = seed_orbits(w, ctx, E, l, Hp, Hq, Hm, Hkeyidx)
            trec.setdefault("alive", []).append(na); trec.setdefault("lam", []).append(nlam)
            for s in reps:
                if len(G.MATS) > MATCAP // 2:
                    # memory: the Game's interned matrices accumulate over components; rebuild it (same prefix)
                    w = W4(1); G = w.G
                    gens = build(w, rec["prog"])
                    ctx = Ctx(w, gens)
                r = component(ctx, E, (l, G.intern(G.snap(SK[s]))), caprows, capcols)
                r["l"] = l; r["seed"] = s
                trec["orbits"].append(r)
                if r["verdict"] == "nonzero":
                    allzero = False
                elif r["verdict"] == "undecided":
                    undecided = True
        out["tests"].append(trec)
    # private-end status per array in the old sense (free_last_test), derived when everything is zero
    ecand = [len(G.E_candidates(gens, l)) for l in range(2)]
    out["Ecand"] = ecand
    out["allzero"] = allzero and not undecided
    out["undecided"] = undecided
    out["nonzero"] = not allzero
    out["sec"] = round(time.time() - t0, 1)
    if out["allzero"]:
        out["PRIV"] = ["none" if ecand[l] == 0 else "no" for l in range(2)]
    return out

ONLY_R = os.environ.get('O3F_GROUP', 'P') == 'R'    # O3F_GROUP=R: o3-build's 384 address maps only

def hdata():
    Hp, Hq, HisR = H_group()
    if ONLY_R:
        Hp, Hq, HisR = Hp[HisR], Hq[HisR], HisR[HisR]
    n = len(Hp)
    Hm = np.zeros((n, 16, 16), dtype=complex)
    Hm[np.arange(n)[:, None], Hp, np.arange(16)[None, :]] = Z4c[Hq]
    Hkeyidx = {(Hp[k].tobytes(), Hq[k].tobytes()): k for k in range(n)}
    return Hp, Hq, HisR, Hm, Hkeyidx

if __name__ == "__main__":
    inp, outp = sys.argv[1], sys.argv[2]
    first = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    step = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    caprows = int(sys.argv[5]) if len(sys.argv) > 5 else 60000
    H = hdata()
    recs = [json.loads(l) for l in open(inp)]
    fo = open(outp, "a")
    for k in range(first, len(recs), step):
        r = test_prefix(recs[k], caprows, caprows, H)
        r["k"] = k
        fo.write(json.dumps(r) + "\n"); fo.flush()
        print(k, r["run"], r["it"], "nE", r["nE"], "Rcls", r["nRclass"], "Pcls", r["nPclass"], "allzero", r["allzero"],
              "undecided", r["undecided"], "nonzero", r["nonzero"], r["sec"], "s", flush=True)
