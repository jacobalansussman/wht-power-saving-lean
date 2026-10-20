"""poolsuff.py -- o3-finish2: the P-test (bases over P = S4 x| Pauli4) against EVERY end move, for kframe prefixes,
by proving that end moves outside the pool and outside a short threat list cannot help (PROOF.md, Lemma P).

Lemma P.  In the P-system for a paid element E (seeds s K, s in P; every other non-free, non-E-free row an equation),
a seed is touched by span columns b u_t (u_t with a K-frame term c K, b in P c^-1) and, if E is in a pool class, by
E-columns.  Let R0_t be a deadness certificate for 1 . u_t in the span-only system (kfD's 'pure' test; sinks there
are all K-frame rows, so R0_t holds no K-frame row and its equations are P-system equations too).  For E outside the
pool, the certificate b R0_t of b u_t (b in P c^-1) fails only if an E-column or an E-free row meets b R0_t, i.e.
E beta g' = b r (r in R0_t, g' a generator term or 1), i.e. E in b r g'^-1 B', whose P-class is that of c^-1 r g'^-1.
So every E outside  pool + T,  T = { P-classes of c^-1 r g'^-1 : c K a K-frame term of a kframe u_t, r in R0_t,
r g'^-1 one move },  has every seed-touching column dead: M_E(P) = 0.  The classes in T are tested directly.
usage: python3 poolsuff.py PREFIXES.jsonl OUT.jsonl [only RUN:IT,...]"""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault('O3F_PRIMES', '1048573,2097133')
from kfD import *
import fulltest as FT
import kfD as _KFD
if os.environ.get('PS_EXACT') == '1':
    _orig_solve = FT.solve_dimM
    def _exact_solve(nrows, ent, tset, primes=None):
        if len(tset) == 1:
            return _KFD.exact_dimT(nrows, ent, tset)
        return _orig_solve(nrows, ent, tset, primes)     # only the informational dim of all seeds stays mod p
    FT.solve_dimM = _exact_solve


def pclass_key(E, Hp, Hq, HisR, Hm):
    Hinv = np.conj(np.transpose(Hm, (0, 2, 1)))
    X = np.matmul(np.conj(E).T[None], Hinv)
    return min(FT.vkeys(X))


def test(rec, Hdata, caprows=40000):
    t0 = time.time()
    Hp, Hq, HisR, Hm, Hkeyidx = Hdata
    S = Sys(rec); G = S.G
    out = {"run": rec.get("run", "").replace(".jsonl", ""), "it": rec.get("it")}
    kgen = [t for t in range(2, len(S.gens)) if any(G.frm(gid) == S.Kfrm for (l, gid) in S.gens[t])]
    # 1. pure certificates
    R0 = {}
    for t in kgen:
        r = dead_test(S, ('S', t, np.eye(16, dtype=complex)), False, caprows)
        if r[0] != "dead":
            out.update(verdict="no-certificate", t=t, status=r[0], sec=round(time.time() - t0, 1))
            return out
        R0[t] = [(l, G.MATS[gid].copy()) for (l, gid) in r[3]]
    out["R0"] = {str(t): len(v) for t, v in R0.items()}
    # 2. threat list
    w = FT.W4(1); G2 = w.G
    gens = FT.build(w, rec["prog"])
    ctx = FT.Ctx(w, gens)
    pool = FT.e_classes(w, gens, Hp, Hq, HisR)
    poolkeys = {kP for (_, _, kP) in pool}
    threats = {}
    for t in kgen:
        cterms = [G2.MATS[gid] for (l, gid) in gens[t] if G2.frm(gid) == ctx.Kfrm]
        for cK in cterms:
            ci = np.conj(G2.mul(cK, np.conj(G2.K).T)).T          # c^-1
            for (l, Mr) in R0[t]:
                for s2 in range(len(gens)):
                    for (l2, g2) in gens[s2]:
                        if l2 != l:
                            continue
                        h = G2.mul(Mr, np.conj(G2.MATS[g2]).T)
                        if decomp1(h) is None:
                            continue
                        E = G2.mul(ci, h)
                        k = pclass_key(E, Hp, Hq, HisR, Hm)
                        if k not in poolkeys:
                            threats.setdefault(k, E)
    if os.environ.get('PS_POOL') == '1':      # also re-test the pool classes themselves (second prime / exact)
        for (E, kR, kP) in pool:
            threats.setdefault(kP, E)
    out["n_pool"] = len(poolkeys); out["n_threat"] = len(threats)
    # 3. test each threat class with the P-test's own machinery (fulltest.component on every seed orbit)
    allzero = True; undecided = False; recs = []
    for k, E in threats.items():
        for l in range(2):
            na, nlam, reps, SK = FT.seed_orbits(w, ctx, E, l, Hp, Hq, Hm, Hkeyidx)
            for s in reps:
                if len(G2.MATS) > FT.MATCAP // 2:
                    w = FT.W4(1); G2 = w.G; gens = FT.build(w, rec["prog"]); ctx = FT.Ctx(w, gens)
                r = FT.component(ctx, E, (l, G2.intern(G2.snap(SK[s]))), caprows, caprows)
                recs.append([l, int(s), r["verdict"], r["status"], r["rows"]])
                if r["verdict"] == "nonzero":
                    allzero = False
                elif r["verdict"] == "undecided":
                    undecided = True
    out["threat_tests"] = recs[:50]
    out["verdict"] = "pool-suffices" if allzero and not undecided else ("NONZERO" if not allzero else "undecided")
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
        r["exact"] = os.environ.get('PS_EXACT') == '1'
        fo.write(json.dumps(r, default=str) + "\n"); fo.flush()
        print(r["run"], r["it"], r["verdict"], "R0", r.get("R0"), "pool", r.get("n_pool"), "threats", r.get("n_threat"),
              r["sec"], "s", flush=True)
