"""fulltestG.py -- o3-finish: the widest form of the private-end test: ANY invertible free mixing of the two end moves.
A 7-move word ending  K x = Q (E6 c6, E7 c7)^T + sigma  with Q in GL_2(Q[B']) (B' = every free letter of the game,
not only those normalising K) exists iff there are rows (al1, al2) with  al1 K x1 + al2 K x2 in span + E_j span
(j = 6, 7) forming an invertible matrix.  So the seeds are ALL rows in K's frame, (l, b K) with b in B', and
    M_E(B') = { (al1, al2) in Q[B']^2 : al1 K x1 + al2 K x2 in span + E span }.
This contains the bases of fulltest.py (P) and fulltestN.py (N): F in GL_2(Q[N]) is Q = K F^-1 K^-1.
Two simplifications make it cheap:
  * every paid element of the game is one move, E = b D1 b', and M_{b D1 b'} = b M_{D1} (the closure system of
    b D1 b' is the image of D1's under left multiplication by b), so ONE test E = D1 covers every E (plus E = None);
  * the seeds touched by D1-columns through a dist-1 term g2 of array l form one orbit of Lambda = B' cap D1 B' D1^-1,
    with representative b3^-1 K where K g2^-1 = b3 D1 b4 (decomposition by bm's halving matching); seeds touched by
    span columns (terms in K's frame) with E = D1 are not enumerated ('kframe': undecided here, decided in P).
usage: python3 -I fulltestG.py PREFIXES.jsonl OUT.jsonl [first] [step] [caprows]"""
import sys, os, json, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fulltest as F
from fulltest import *

class CtxG(F.Ctx):
    def seed_of(self, gid):
        return self.G.frm(gid) == self.Kfrm

def decomp(G, h):
    """h = b3 D1 b4 with b3, b4 free, or None"""
    Ec = G.compat(h)
    rng = random.Random(1)
    D1i = np.conj(G.D1).T
    for M in G.rand_matchings(Ec, rng, 30, 200):
        m = G.back_step(h, M, Ec)
        if G.is_monomial(m):
            b3 = G.mul(h, np.conj(m).T, D1i)
            if G.is_monomial(b3) and np.allclose(G.mul(b3, G.D1, m), h):
                return b3, m
    return None

def test_prefixG(rec, caprows):
    w = W4(1); G = w.G
    gens = build(w, rec["prog"])
    ctx = CtxG(w, gens)
    t0 = time.time()
    out = {"run": rec.get("run"), "it": rec.get("it"), "tests": []}
    allzero = True; undecided = False; kframe = False
    d1 = {0: [], 1: []}
    for g in gens:
        for (l, gid) in g:
            if w.dist(gid) == 1 and gid not in d1[l]:
                d1[l].append(gid)
    nonone = 0
    for (nm, E) in (("None", None), ("D1", G.D1)):
        trec = {"E": nm, "orbits": []}
        for l in range(2):
            reps = []
            if E is None:
                if (ctx.Kfrm, l) in ctx.index:
                    reps.append(G.KID)
            else:
                if (ctx.Kfrm, l) in ctx.index:
                    kframe = True; trec["kframe"] = True
                for gid2 in d1[l]:
                    h = G.mul(G.K, np.conj(G.MATS[gid2]).T)
                    if not (G.rowsupp(h) == 2).all():
                        continue
                    d = decomp(G, h)
                    if d is None:
                        nonone += 1; continue
                    b3, b4 = d
                    reps.append(G.intern(G.mul(np.conj(b3).T, G.K)))
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
    out["d1terms"] = [len(d1[0]), len(d1[1])]; out["not_one_move"] = nonone
    out["allzero"] = allzero and not undecided and not kframe
    out["undecided"] = undecided; out["kframe"] = kframe; out["nonzero"] = not allzero
    out["sec"] = round(time.time() - t0, 1)
    return out

if __name__ == "__main__":
    inp, outp = sys.argv[1], sys.argv[2]
    first = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    step = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    caprows = int(sys.argv[5]) if len(sys.argv) > 5 else 40000
    recs = [json.loads(l) for l in open(inp)]
    fo = open(outp, "a")
    for k in range(first, len(recs), step):
        r = test_prefixG(recs[k], caprows)
        r["k"] = k
        fo.write(json.dumps(r) + "\n"); fo.flush()
        print(k, r["run"], r["it"], "allzero", r["allzero"], "undecided", r["undecided"], "kframe", r["kframe"],
              "nonzero", r["nonzero"], "d1", r["d1terms"], r["sec"], "s", flush=True)
