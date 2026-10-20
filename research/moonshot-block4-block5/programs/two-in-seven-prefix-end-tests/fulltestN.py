"""fulltestN.py -- o3-finish: the any-basis test of fulltest.py with the FULL normaliser
    N = { h free : K h K^-1 free }  (20,643,840 letters: affine address maps with the twist i^(wt(pi y) - wt y),
                                       Z_b, the scalar letter iI; it contains P = S4 x| Pauli4 and R)
as basis group: F in GL_2(Q[N]) is the largest group of free basis changes for which K F K^-1 is again free, i.e.
for which Lemma 1 + Cor. 1b turn the end into two private ends.

Seeds (l, K h), h in N, are never enumerated.  A seed is touched by an E-column only through a generator term g2 of
array l with E^-1 K h = g_t h in the frame of g2 (g_t = E^-1 K, itself a dist-1 term), and the seeds touched through
one g2 form ONE orbit of Lambda_E = N cap E B' E^-1 (if g_t h0 and g_t h1 are both in B' g2 then h1 h0^-1 lies in
N cap g_t^-1 B' g_t = K^-1 Lambda_E K).  So one closure per (class of E, term g2) with a representative h0 found by
matching supports under AGL(4,2) (322,560 maps) and then the 64 phases i^c Z_b.  Seeds touched by span columns (a term
of array l in K's frame) are all of N: for E = None that is one orbit (Lambda = N); for E != None this case is not
enumerated and is reported as 'kframe' (undecided in N, still decided in P by fulltest.py).
usage: python3 -I fulltestN.py PREFIXES.jsonl OUT.jsonl [first] [step] [caprows]"""
import sys, os, json, time, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fulltest as F
from fulltest import *

def agl():
    mats = []
    for cols in itertools.product(range(1, 16), repeat=4):
        # invertible iff the 4 columns span F_2^4
        span = {0}
        for c in cols:
            span |= {s ^ c for s in span}
        if len(span) == 16:
            mats.append(cols)
    out = []
    for cols in mats:
        lin = [0] * 16
        for y in range(16):
            v = 0
            for j in range(4):
                if (y >> j) & 1: v ^= cols[j]
            lin[y] = v
        for a in range(16):
            out.append([lin[y] ^ a for y in range(16)])
    return np.array(out, dtype=np.int64)

WT = np.array([bin(y).count('1') for y in range(16)])
PI = None

def supp_masks(M):
    return (np.abs(M) > 1e-9)

def canon(S):
    """AGL-canonical form of the row-support family S (16 x 16 bool): lexmin over pi of the sorted row masks of S[:, pi]"""
    wts = (1 << np.arange(16)).astype(np.int32)
    m = np.empty((PI.shape[0], 16), dtype=np.int32)
    for r in range(16):
        m[:, r] = (S[r][PI].astype(np.int32) * wts[None, :]).sum(axis=1)
    m.sort(axis=1)
    o = np.lexsort(m.T[::-1])
    return m[o[0]].tobytes(), m

def in_N(G, m):
    return G.is_monomial(G.mul(G.K, m, np.conj(G.K).T))

def find_h(G, gt, g2, mt, sorted2):
    """h in N with frame(gt h) == frame(g2), or None.  mt: (n, 16) sorted masks of gt's supports under every pi"""
    cand = np.nonzero((mt == sorted2[None, :]).all(axis=1))[0]
    k2 = F.vkeys(g2[None])[0]
    for ci in cand:
        p = PI[ci]
        tw = (WT[p] - WT) % 4
        hs = []
        for c in range(4):
            for b in range(16):
                q = (tw + c + 2 * np.array([dot(b, int(p[y])) for y in range(16)])) % 4
                hs.append(F.pq_mat(p, q))
        Hs = np.array(hs)
        ks = F.vkeys(np.matmul(gt[None], Hs))
        for j, kk in enumerate(ks):
            if kk == k2:
                h = G.snap(Hs[j])
                if in_N(G, h):
                    return h
    return None

class CtxN(F.Ctx):
    def seed_of(self, gid):
        v = self._seedmemo.get(gid)
        if v is None:
            G = self.G
            v = False
            if G.frm(gid) == self.Kfrm:
                m = G.mul(G.MATS[gid], self.Ki)
                v = G.is_monomial(m) and in_N(G, m)
            self._seedmemo[gid] = v
        return v

def test_prefixN(rec, caprows, Hdata):
    w = W4(1); G = w.G
    gens = build(w, rec["prog"])
    ctx = CtxN(w, gens)
    Hp, Hq, HisR, Hm, Hkeyidx = Hdata
    t0 = time.time()
    ecl = F.e_classes(w, gens, Hp, Hq, HisR)
    Pcls = {}
    for (E, kR, kP) in ecl:
        Pcls.setdefault(kP, E)
    tests = [("None", None)] + [("E%d" % k, E) for k, E in enumerate(Pcls.values())]
    # dist-1 terms per array, with support canonical forms
    d1 = {0: [], 1: []}
    cache = {}
    for g in gens:
        for (l, gid) in g:
            if w.dist(gid) == 1 and gid not in [x for x in d1[l]]:
                d1[l].append(gid)
    def info(gid):
        if gid not in cache:
            S = supp_masks(G.MATS[gid])
            c, m = canon(S)
            s2 = np.sort((S.astype(np.int32) * (1 << np.arange(16)).astype(np.int32)[None, :]).sum(axis=1))
            cache[gid] = (c, m, s2)
        return cache[gid]
    out = {"run": rec.get("run"), "it": rec.get("it"), "nPclass": len(Pcls), "tests": []}
    allzero = True; undecided = False; kframe = False
    for (nm, E) in tests:
        trec = {"E": nm, "orbits": []}
        for l in range(2):
            reps = []
            if E is None:
                if (ctx.Kfrm, l) in ctx.index:
                    reps.append(G.KID)
            else:
                if (ctx.Kfrm, l) in ctx.index:
                    kframe = True; trec["kframe"] = True
                gt = G.mul(np.conj(E).T, G.K)
                gtid = G.intern(gt)
                ct, mt, _ = info(gtid)
                for gid2 in d1[l]:
                    c2, _, s2 = info(gid2)
                    if c2 != ct:
                        continue
                    h = find_h(G, gt, G.MATS[gid2], mt, s2)
                    if h is not None:
                        reps.append(G.intern(G.mul(G.K, h)))
            reps = list(dict.fromkeys(reps))
            for sid in reps:
                r = F.component(ctx, E, (l, sid), caprows, caprows)
                r["l"] = l
                if r["verdict"] == "nonzero":
                    allzero = False
                elif r["verdict"] == "undecided":
                    undecided = True
                trec["orbits"].append(r)
        out["tests"].append(trec)
    out["allzero"] = allzero and not undecided and not kframe
    out["undecided"] = undecided; out["kframe"] = kframe; out["nonzero"] = not allzero
    out["sec"] = round(time.time() - t0, 1)
    return out

if __name__ == "__main__":
    inp, outp = sys.argv[1], sys.argv[2]
    first = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    step = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    caprows = int(sys.argv[5]) if len(sys.argv) > 5 else 40000
    PI = agl(); F.PI = PI
    globals()['PI'] = PI
    H = F.hdata()
    recs = [json.loads(l) for l in open(inp)]
    fo = open(outp, "a")
    for k in range(first, len(recs), step):
        r = test_prefixN(recs[k], caprows, H)
        r["k"] = k
        fo.write(json.dumps(r) + "\n"); fo.flush()
        print(k, r["run"], r["it"], "Pcls", r["nPclass"], "allzero", r["allzero"], "undecided", r["undecided"],
              "kframe", r["kframe"], "nonzero", r["nonzero"], r["sec"], "s", flush=True)
