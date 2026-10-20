"""kfG.py -- o3-finish2 (14 Oct 2026): the WIDEST end test (o3-finish's G-test) extended to kframe prefixes.

Question for a 5-move prefix (generators y1, y2, u1..u5): is there an end
    K x = Q (E6 c6, E7 c7)^T + s,   Q in GL_2(Q[B']) (B' = every free letter), E6, E7 ANY paid elements, c6, c7, s in S5?
o3-finish (proved, skeptic-upheld): every paid element is b D1 b', M_{b D1 b'} = b M_D1, and a solution needs
dim eps(M_D1) = 2 (eps = coefficient sum, a ring map; eps(b v) = eps(v)).  Here M_D1 = { (a1, a2) in Q[B']^2 :
a1 K x1 + a2 K x2 in S5 + D1 S5 }, seeds = ALL rows of K's frame.

Seeds as sinks.  A seed row (l, bK) carries its own free unknown t, so it imposes no equation: the system is
'every NON-seed row cancels', and eps(M_D1) is the image of its solution space under
    o_l(x) = sum over seed rows of input l of (row . x)        (l = 1, 2).
Components are therefore taken through non-seed rows only, and
    eps(M_D1) = sum over components C touching a seed of o(ker A_C).
Free rows (monomial) and E-free rows (D1^-1 g monomial) are covered by single-entry columns (raw input; D1 b y_l) and are
dropped exactly.

Which components (the symmetry): Lambda = B' cap D1 B' D1^-1 acts on the left and permutes components (lam D1 beta =
D1 (D1^-1 lam D1) beta); B' itself permutes the components of the span-only system.  A component touching a seed is
  (P) PURE: span columns only, no E-free row: a B'-translate b.N_t of the span-only component N_t of the column u_t
      (u_t a generator with a term in K's frame); same eps image as N_t;
  (S) SPECIAL: it holds a D1 column (or an E-free row).  Either a D1 column touches a seed directly -- then it is a
      Lambda-translate of the component of D1 b4 u_s, K g^-1 = b3 D1 b4 for a term g of u_s ('direct starts') -- or the
      first D1 column on a path from a seed-touching span column b u_t touches a non-seed row b r with r in N_t; then
      D1 beta g' = b r for a term g' of some generator (g' = 1 for a raw input), i.e. r g'^-1 = c1 D1 c2 one move and
      b in Lambda c1^-1: a Lambda-translate of the component of the span column c1^-1 u_t ('translate starts').
So the starts are finite once every N_t is complete.  A component explored only up to a cap is replaced by its
relaxation (processed rows; every column touching them): fewer equations, larger image -> sound for 'dim <= 1'.
If some N_t is incomplete (capped) the translate starts are incomplete: the verdict is then UNDECIDED.

Linear algebra mod two primes (screening); 'excluded' needs dim <= 1 mod both, and the images of all components are
compared (sum of images).  Exact confirmation: kfG_exact.py.
usage: python3 kfG.py PREFIXES.jsonl OUT.jsonl [first] [step] [caprows] [only RUN:IT,...]"""
import sys, os, json, time, collections
PYLIB = ''
if os.path.isdir(PYLIB):
    sys.path.insert(0, PYLIB)
HERE = os.path.dirname(os.path.abspath(__file__))
V3 = os.path.join(HERE, 'v3')
if not os.path.isdir(V3):
    raise SystemExit('v3/ not found next to kfG.py')
sys.path.insert(0, V3)
import numpy as np
from fractions import Fraction
from o3lib import W4, build
import bm as BM

PRIMES = (1048573, 2097133)
A_ = (1 + 1j) / 2; B_ = (1 - 1j) / 2
MATCAP = int(os.environ.get('KFG_MATCAP', 400000))
VERBOSE = os.environ.get('KFG_VERBOSE', '0') == '1'


def rankp(A, pr):
    if A.size == 0:
        return 0
    if A.shape[1] > A.shape[0]:
        A = A.T
    A = np.ascontiguousarray(np.mod(A, pr).astype(np.uint64))
    R, C = A.shape
    piv = np.zeros(max(1, min(R, C)), dtype=np.int64)
    return int(BM._lib.rref(A.ctypes.data, R, C, pr, piv.ctypes.data, 0))


def decomp1(h):
    """h = c1 D1 c2 with c1, c2 Z4-monomials, or None (deterministic; see NOTE: rows of support 2, paired by support,
    partner ratio = - ratio)"""
    n = h.shape[0]
    nz = np.abs(h) > 1e-9
    if not (nz.sum(axis=1) == 2).all():
        return None
    sup = {}
    for x in range(n):
        sup.setdefault(tuple(np.nonzero(nz[x])[0]), []).append(x)
    if any(len(v) != 2 for v in sup.values()):
        return None
    cols = sorted(c for k in sup for c in k)
    if cols != list(range(n)):
        return None
    c1 = np.zeros((n, n), dtype=complex); c2 = np.zeros((n, n), dtype=complex)
    for k, ((ca, cb), (x, x2)) in enumerate(sup.items()):
        rho = h[x, cb] / h[x, ca]
        p = 1j * rho
        c2[2 * k, ca] = 1; c2[2 * k + 1, cb] = p
        lam = h[x, ca] / A_; lam2 = h[x2, ca] / B_
        c1[x, 2 * k] = lam; c1[x2, 2 * k + 1] = lam2
    if not (BM.Game.is_monomial(c1) and BM.Game.is_monomial(c2)):
        return None
    if not np.allclose(c1 @ (A_ * np.eye(n) + B_ * np.eye(n)[[y ^ 1 for y in range(n)]]) @ c2, h, atol=1e-9):
        return None
    return c1, c2


class Sys:
    def __init__(self, rec):
        self.w = W4(1); G = self.G = self.w.G
        self.gens = build(self.w, rec["prog"])
        self.index = {}
        for i, g in enumerate(self.gens):
            for (ll, gid) in g:
                self.index.setdefault((G.frm(gid), ll), []).append((i, gid))
        self.Kfrm = G.frm(G.KID)
        self.D1 = G.D1; self.D1i = np.conj(G.D1).T
        self._ef = {}

    def seed(self, gid):
        return self.G.frm(gid) == self.Kfrm

    def efree(self, gid):
        v = self._ef.get(gid)
        if v is None:
            v = self.G.is_monomial(self.G.mul(self.D1i, self.G.MATS[gid])); self._ef[gid] = v
        return v

    def col_entries(self, kind, i, bmat, withE):
        """rows of the column b gen_i ('S') or D1 b gen_i ('E'): list of ((l, gid) or ('seed', l), coef)"""
        G = self.G
        P2 = bmat if kind == 'S' else G.mul(self.D1, bmat)
        out = []
        for (l2, g2), co in self.gens[i].items():
            g3 = G.intern(G.mul(P2, G.MATS[g2]))
            if G.is_free(g3):
                continue
            if withE and self.efree(g3):
                continue
            if self.seed(g3):
                out.append((('seed', l2), co))
            else:
                out.append(((l2, g3), co))
        return out

    def component(self, starts, withE, caprows):
        """BFS through non-seed rows from the given columns [(kind, i, bmat)]; returns dict with rows, cols, entries"""
        G = self.G
        cols = {}; rows = {}; ent = []        # ent: (rowkey, colidx, coef); rowkey ('seed', l) for the outputs
        q = collections.deque(); processed = []
        def addcol(kind, i, bmat):
            b = G.intern(bmat)
            ck = (kind, i, b)
            if ck in cols:
                return
            c = len(cols); cols[ck] = c
            for (rk, co) in self.col_entries(kind, i, bmat, withE):
                ent.append((rk, c, co))
                if rk[0] != 'seed' and rk not in rows:
                    rows[rk] = len(rows); q.append(rk)
        for (kind, i, bmat) in starts:
            addcol(kind, i, bmat)
        capped = False
        while q:
            if len(processed) >= caprows or len(G.MATS) > MATCAP:
                capped = True; break
            rk = q.popleft(); processed.append(rk)
            (l, gid) = rk
            Mg = G.MATS[gid]
            for (i, gterm) in self.index.get((G.frm(gid), l), ()):
                if i < 2:
                    continue
                addcol('S', i, G.mul(Mg, np.conj(G.MATS[gterm]).T))
            if withE:
                hm = G.mul(self.D1i, Mg)
                for (i, gterm) in self.index.get((G.frame_key(hm), l), ()):
                    if i < 2:
                        continue
                    addcol('E', i, G.mul(hm, np.conj(G.MATS[gterm]).T))
        return {"cols": cols, "rows": rows, "ent": ent, "processed": processed, "capped": capped}

    def image(self, comp, pr, maxcore=1.5e7):
        """basis of o(ker A) mod pr (list of 2-vectors); None if the component is capped (unprocessed rows: the
        relaxation would leave unseen columns free) or its peeled core is too large"""
        if comp["capped"]:
            return None
        rows = collections.defaultdict(dict); colsd = collections.defaultdict(set); O = collections.defaultdict(dict)
        for (rk, c, co) in comp["ent"]:
            v = (co.numerator % pr) * pow(co.denominator % pr, pr - 2, pr) % pr
            if rk[0] == 'seed':
                O[c][rk[1]] = (O[c].get(rk[1], 0) + v) % pr
            else:
                nv = (rows[rk].get(c, 0) + v) % pr
                if nv:
                    rows[rk][c] = nv; colsd[c].add(rk)
                else:
                    rows[rk].pop(c, None); colsd[c].discard(rk)
        outc = {c for c, d in O.items() if any(d.values())}
        changed = True
        while changed:
            changed = False
            for r in list(rows.keys()):
                rr = rows.get(r)
                if rr is None:
                    continue
                if not rr:
                    del rows[r]; changed = True; continue
                if len(rr) == 1:
                    (c,) = rr.keys()
                    for r2 in list(colsd[c]):
                        rows[r2].pop(c, None)
                    colsd[c] = set(); outc.discard(c)
                    del rows[r]; changed = True
            for c in list(colsd.keys()):
                if c in outc or len(colsd[c]) != 1:
                    continue
                (r,) = colsd[c]
                if r in rows:
                    for c2 in list(rows[r].keys()):
                        colsd[c2].discard(r)
                    del rows[r]
                colsd[c] = set(); changed = True
        vecs = []
        freeout = [c for c in outc if not colsd.get(c)]
        for c in freeout:
            vecs.append((O[c].get(0, 0), O[c].get(1, 0)))
        rl = list(rows.keys())
        cl = sorted({c for r in rl for c in rows[r]})
        if not rl or not any(c in outc for c in cl):
            return vecs
        if len(rl) * len(cl) > maxcore:
            return None
        ci = {c: k for k, c in enumerate(cl)}
        A = np.zeros((len(rl), len(cl)), dtype=np.float64)
        for i, r in enumerate(rl):
            for c, v in rows[r].items():
                A[i, ci[c]] = v
        Om = np.zeros((2, len(cl)), dtype=np.float64)
        for c in cl:
            if c in outc:
                Om[0, ci[c]] = O[c].get(0, 0); Om[1, ci[c]] = O[c].get(1, 0)
        rA = rankp(A, pr); rAO = rankp(np.vstack([A, Om]), pr)
        d = rAO - rA
        if d == 2:
            return vecs + [(1, 0), (0, 1)]
        if d == 1:
            Nm = BM_nullspace(A, pr)
            img = np.mod(Om @ Nm, pr) if Nm.size else np.zeros((2, 0))
            for j in range(img.shape[1]):
                v = (int(img[0, j]), int(img[1, j]))
                if v != (0, 0):
                    vecs.append(v); break
        return vecs


def BM_nullspace(A, pr):
    R, C = A.shape
    A2 = np.ascontiguousarray(np.mod(A, pr).astype(np.uint64))
    piv = np.zeros(max(1, min(R, C)), dtype=np.int64)
    r = int(BM._lib.rref(A2.ctypes.data, R, C, pr, piv.ctypes.data, 1))
    Rr = A2[:r].astype(np.int64); pv = list(piv[:r])
    free = [c for c in range(C) if c not in set(pv)]
    Nm = np.zeros((C, len(free)), dtype=np.float64)
    for t, fc in enumerate(free):
        Nm[fc, t] = 1
        for i2, pc in enumerate(pv):
            Nm[pc, t] = (-Rr[i2, fc]) % pr
    return Nm


def span_dim(vecs, pr):
    vs = [v for v in vecs if v != (0, 0)]
    if not vs:
        return 0
    a = vs[0]
    for b in vs[1:]:
        if (a[0] * b[1] - a[1] * b[0]) % pr:
            return 2
    return 1


def test(rec, caprows=40000):
    t0 = time.time()
    S = Sys(rec); G = S.G
    gens = S.gens
    out = {"run": rec.get("run", "").replace(".jsonl", ""), "it": rec.get("it")}
    kgen = sorted({i for i, g in enumerate(gens) if i >= 2 and any(G.frm(gid) == S.Kfrm for (l, gid) in g)})
    out["kgens"] = kgen
    comps = []          # (label, comp)
    # (P) span-only components of the kframe generators
    Nt = {}
    for t in kgen:
        c = S.component([('S', t, np.eye(16, dtype=complex))], False, caprows)
        Nt[t] = c; comps.append((('N', t), c))
    if VERBOSE:
        print("   N_t", {t: (len(Nt[t]["rows"]), Nt[t]["capped"]) for t in kgen}, round(time.time() - t0, 1),
              file=sys.stderr, flush=True)
    out["N_rows"] = {str(t): len(Nt[t]["rows"]) for t in kgen}
    out["N_capped"] = any(Nt[t]["capped"] for t in kgen)
    # (S) starts
    starts = {}
    # direct: D1 b4 u_s with D1 b4 g in K's frame
    for s in range(2, len(gens)):
        for (l, gid) in gens[s]:
            h = G.mul(G.K, np.conj(G.MATS[gid]).T)
            d = decomp1(h)
            if d is None:
                continue
            b3, b4 = d
            key = ('dir', s, G.frame_key(G.mul(G.D1, b4)))
            starts.setdefault(key, ('E', s, b4))
    # translates: rows r of N_t, terms g' of any generator (raw inputs too) on the same input, r g'^-1 one move
    nterm = 0
    for t in kgen:
        for (l, gid) in Nt[t]["rows"]:
            Mr = G.MATS[gid]
            for s in range(len(gens)):
                for (l2, g2) in gens[s]:
                    if l2 != l:
                        continue
                    h = G.mul(Mr, np.conj(G.MATS[g2]).T)
                    d = decomp1(h)
                    if d is None:
                        continue
                    c1, c2 = d
                    nterm += 1
                    bt = np.conj(c1).T
                    key = ('tr', t, G.frame_key(G.mul(S.D1i, bt)))
                    starts.setdefault(key, ('S', t, bt))
    out["n_starts"] = len(starts); out["n_pairs"] = nterm
    seen_lam = set()        # Lambda-orbit keys of every column of the special components explored so far
    def lamkey(kind, i, bmat):
        if kind == 'S':
            return ('tr', i, G.frame_key(G.mul(S.D1i, bmat)))
        return ('dir', i, G.frame_key(G.mul(G.D1, bmat)))
    nspec = 0; spec_capped = 0; spec_rows = []; nskip = 0
    for key, st in starts.items():
        if key in seen_lam:
            nskip += 1
            continue
        c = S.component([st], True, caprows)
        for (kind, i, b) in c["cols"].keys():
            seen_lam.add(lamkey(kind, i, G.MATS[b]))
        comps.append((key, c)); nspec += 1
        spec_capped += c["capped"]; spec_rows.append(len(c["rows"]))
        if len(G.MATS) > MATCAP:
            break
        if VERBOSE:
            print("   special", nspec, "of", len(starts), "rows", len(c["rows"]), "cols", len(c["cols"]), "capped",
                  c["capped"], "mats", len(G.MATS), round(time.time() - t0, 1), file=sys.stderr, flush=True)
    out["n_special"] = nspec; out["special_capped"] = spec_capped; out["n_lam_skipped"] = nskip
    out["special_rows_max"] = max(spec_rows, default=0)
    out["incomplete_starts"] = len(G.MATS) > MATCAP and nspec < len(starts)
    res = {}
    for pr in PRIMES:
        vecs = []; big = 0
        for (lab, c) in comps:
            im = S.image(c, pr)
            if im is None:
                big += 1
                continue
            vecs += im
        res[pr] = (span_dim(vecs, pr), big, [v for v in vecs if v != (0, 0)][:6])
    out["dim_eps"] = [res[p][0] for p in PRIMES]
    out["big"] = [res[p][1] for p in PRIMES]
    out["vecs"] = [[list(map(int, v)) for v in res[p][2]] for p in PRIMES]
    complete = not out["N_capped"] and not out["incomplete_starts"] and max(out["big"]) == 0
    out["complete_starts"] = complete
    dmax = max(out["dim_eps"])
    out["verdict"] = ("excluded" if dmax <= 1 and complete else
                      ("undecided" if dmax <= 1 else "eps-dim-2"))
    out["sec"] = round(time.time() - t0, 1)
    return out


if __name__ == "__main__":
    inp, outp = sys.argv[1], sys.argv[2]
    first = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    step = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    caprows = int(sys.argv[5]) if len(sys.argv) > 5 else 40000
    only = None
    if len(sys.argv) > 6 and sys.argv[6]:
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
        print(kk, r["run"], r["it"], "kgens", r["kgens"], "N", r["N_rows"], "Ncap", r["N_capped"], "starts",
              r["n_starts"], "special", r["n_special"], "cap", r["special_capped"], "max", r["special_rows_max"],
              "dim", r["dim_eps"], r["verdict"], r["sec"], "s", flush=True)
