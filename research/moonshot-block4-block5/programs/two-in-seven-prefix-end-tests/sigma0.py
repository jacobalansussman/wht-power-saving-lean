"""sigma0.py -- o3-finish2: Criterion 0 of the augmentation lemma (exact, over Q, no closure).
Every column of the widest end system that touches a seed (a row in K's frame) is either a span column b u_t (u_t a
generator with a term in K's frame) or a D1 column D1 beta u_s whose beta lies in a coset H_D b4 with K g^-1 = b3 D1 b4
for a term g of u_s.  Its output (K-frame coefficient sums per input) is sigma(u_t), resp. sigma(D1 b4 u_s) (constant on
the coset).  eps(M_D1) lies in the Q-span of these finitely many vectors, so if they span <= 1 dimension the prefix is
excluded for every end move and every mixing over Q[B'] (Lemma A of PROOF.md).
usage: python3 sigma0.py PREFIXES.jsonl OUT.jsonl"""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from kfG import *


def sigma(S, gen):
    G = S.G
    s = [Fraction(0), Fraction(0)]
    for (l, gid), co in gen.items():
        if G.frm(gid) == S.Kfrm:
            s[l] += Fraction(co)
    return s


def rank2(vecs):
    vs = [v for v in vecs if v != [0, 0]]
    if not vs:
        return 0
    a = vs[0]
    return 2 if any(a[0] * b[1] - a[1] * b[0] != 0 for b in vs[1:]) else 1


def crit0(rec):
    S = Sys(rec); G = S.G
    vecs = []; info = []
    for t in range(2, len(S.gens)):
        v = sigma(S, S.gens[t])
        if v != [0, 0]:
            vecs.append(v); info.append(("span", t, [str(x) for x in v]))
    seen = set()
    for s in range(2, len(S.gens)):
        for (l, gid) in S.gens[s]:
            d = decomp1(G.mul(G.K, np.conj(G.MATS[gid]).T))
            if d is None:
                continue
            b4 = d[1]
            key = (s, G.frame_key(G.mul(G.D1, b4)))
            if key in seen:
                continue
            seen.add(key)
            col = G.lmul(G.mul(G.D1, b4), S.gens[s])
            v = sigma(S, col)
            vecs.append(v); info.append(("D1", s, [str(x) for x in v]))
    return rank2(vecs), info


if __name__ == "__main__":
    recs = [json.loads(l) for l in open(sys.argv[1])]
    fo = open(sys.argv[2], "w")
    for k, rec in enumerate(recs):
        t0 = time.time()
        d, info = crit0(rec)
        r = {"run": rec.get("run", "").replace(".jsonl", ""), "it": rec.get("it"), "k": k, "dim0": d,
             "outputs": info, "sec": round(time.time() - t0, 2)}
        fo.write(json.dumps(r) + "\n"); fo.flush()
        print(k, r["run"], r["it"], "dim0", d, "outputs", info[:6], r["sec"], flush=True)
