"""poolsuff2.py -- o3-finish2: Lemma P with the threat classes decided column by column (much cheaper than testing
every seed orbit of a threat class).

For a threat class E (outside the pool), seeds s K (s in P) are touched only by span columns b u_t, b = p c^-1 (p in P,
c K a K-frame term of u_t).  The certificate R0_t transfers to b u_t unless E^-1 b r g''^-1 is free for some r in R0_t
and some generator term g'' on r's input (g'' = 1 allowed).  Those b are found by one batched test over all 24,576 p.
For each such b the column b u_t is tested for deadness directly in the P-system of E (seeds = P K are sinks, other
K-frame rows are equations, E-free rows dropped, E-columns E beta g).  If every such column is dead, every seed
coordinate is a combination of dead columns and M_E(P) = 0.
usage: python3 poolsuff2.py PREFIXES.jsonl OUT.jsonl [only RUN:IT,...]"""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault('O3F_PRIMES', '1048573,2097133')
from kfD import *
import kfD
import fulltest as FT
from poolsuff import pclass_key


class PSys(Sys):
    """the P-system of a paid element E (dead_test works unchanged: D1 -> E)"""
    def __init__(self, rec, E):
        super().__init__(rec)
        self.D1 = E; self.D1i = np.conj(E).T
        self._sd = {}

    def seed(self, gid):
        v = self._sd.get(gid)
        if v is None:
            G = self.G
            v = False
            if G.frm(gid) == self.Kfrm:
                m = G.mul(G.MATS[gid], np.conj(G.K).T)
                if G.is_monomial(m):
                    p, q = FT.mono_pq(G, m)
                    v = FT.in_H(p, q)
            self._sd[gid] = v
        return v


def is_mono_batch(X):
    return FT.is_mono_stack(X)


def test(rec, Hdata, caprows=40000):
    t0 = time.time()
    Hp, Hq, HisR, Hm, Hkeyidx = Hdata
    S = Sys(rec); G = S.G
    out = {"run": rec.get("run", "").replace(".jsonl", ""), "it": rec.get("it")}
    kgen = [t for t in range(2, len(S.gens)) if any(G.frm(gid) == S.Kfrm for (l, gid) in S.gens[t])]
    R0 = {}
    for t in kgen:
        r = dead_test(S, ('S', t, np.eye(16, dtype=complex)), False, caprows)
        if r[0] != "dead":
            out.update(verdict="no-certificate", t=t, status=r[0], sec=round(time.time() - t0, 1))
            return out
        R0[t] = [(l, G.MATS[gid].copy()) for (l, gid) in r[3]]
    out["R0"] = {str(t): len(v) for t, v in R0.items()}
    w = FT.W4(1); G2 = w.G
    gens = FT.build(w, rec["prog"])
    ctx = FT.Ctx(w, gens)
    poolkeys = {kP for (_, _, kP) in FT.e_classes(w, gens, Hp, Hq, HisR)}
    termmats = {l: [G2.MATS[g2] for s2 in range(len(gens)) for (l2, g2) in gens[s2] if l2 == l] for l in range(2)}
    threats = {}
    for t in kgen:
        for cK in [G2.MATS[gid] for (l, gid) in gens[t] if G2.frm(gid) == ctx.Kfrm]:
            ci = np.conj(G2.mul(cK, np.conj(G2.K).T)).T
            for (l, Mr) in R0[t]:
                for Mg in termmats[l]:
                    h = G2.mul(Mr, np.conj(Mg).T)
                    if decomp1(h) is None:
                        continue
                    E = G2.mul(ci, h)
                    k = pclass_key(E, Hp, Hq, HisR, Hm)
                    if k not in poolkeys:
                        threats.setdefault(k, E)
    out["n_pool"] = len(poolkeys); out["n_threat"] = len(threats)
    print("  ", out["run"], out["it"], "R0", out["R0"], "pool", len(poolkeys), "threat classes", len(threats),
          round(time.time() - t0, 1), "s", file=sys.stderr, flush=True)
    # the one-move products r g''^-1 (only these can meet an end move: E beta g'' = b r forces r g''^-1 one move)
    Mlist = {}
    for t in kgen:
        Mlist[t] = []
        for (l, Mr) in R0[t]:
            for Mg in termmats[l] + [np.eye(16, dtype=complex)]:
                M = G2.mul(Mr, np.conj(Mg).T)
                if decomp1(M) is not None:
                    Mlist[t].append((l, M))
    out["n_onemove"] = {str(t): len(v) for t, v in Mlist.items()}
    # per threat class: the columns b u_t whose certificate may fail, tested directly
    report = []; ok = True
    for k, E in threats.items():
        Ei = np.conj(E).T
        cols = []
        for t in kgen:
            for cK in [G2.MATS[gid] for (l, gid) in gens[t] if G2.frm(gid) == ctx.Kfrm]:
                c = G2.mul(cK, np.conj(G2.K).T); ci = np.conj(c).T
                B = np.matmul(Hm, ci[None])                    # b = p c^-1 for all p in P
                bad = np.zeros(len(Hm), dtype=bool)
                for (l, M) in Mlist[t]:
                    if True:
                        X = np.matmul(np.matmul(Ei[None], B), M[None])
                        bad |= is_mono_batch(X)
                for j in np.nonzero(bad)[0]:
                    cols.append((t, G2.snap(B[j])))
        print("   threat class: columns to test", len(cols), round(time.time() - t0, 1), "s", file=sys.stderr, flush=True)
        n_dead = 0; fail = None
        PS = None
        seenb = set()
        for (t, b) in cols:
            kb = (t, G2.key(b))
            if kb in seenb:
                continue
            seenb.add(kb)
            if PS is None or len(PS.G.MATS) > MATCAP // 2:
                PS = PSys(rec, E)
            r = dead_test(PS, ('S', t, b), True, caprows)
            print("     column u%d: %s rows %d cols %d" % (t - 1, r[0], r[1], r[2]), round(time.time() - t0, 1), "s",
                  file=sys.stderr, flush=True)
            if r[0] == "dead":
                n_dead += 1
            else:
                fail = [t, r[0], r[1], r[2]]; ok = False
                break
        report.append({"n_cols": len(seenb), "dead": n_dead, "fail": fail})
        if not ok:
            break
    out["threats"] = report
    out["verdict"] = "pool-suffices" if ok else "undecided"
    out["sec"] = round(time.time() - t0, 1)
    return out


if __name__ == "__main__":
    only = None
    if len(sys.argv) > 3 and sys.argv[3]:
        only = {(x.split(':')[0], int(x.split(':')[1])) for x in sys.argv[3].split(',')}
    H = FT.hdata()
    fo = open(sys.argv[2], "a")
    for l in open(sys.argv[1]):
        rec = json.loads(l)
        key = (rec.get("run", "").replace(".jsonl", ""), rec.get("it"))
        if only is not None and key not in only:
            continue
        r = test(rec, H)
        fo.write(json.dumps(r, default=str) + "\n"); fo.flush()
        print(r["run"], r["it"], r["verdict"], "R0", r.get("R0"), "pool", r.get("n_pool"), "threats", r.get("n_threat"),
              r.get("threats"), r["sec"], "s", flush=True)
