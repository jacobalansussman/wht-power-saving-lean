"""bf7b.py (mw2-b5beyond, resumed 14 Oct ~15:00): special-position butterflies built TRUNK FIRST.
Array 1 gets a full single-letter geodesic chain s1 -> t1 -> t2 -> t3 -> e1 (letters b1, a1, c2, c3, beta1).  Then array 2's
two-letter start q = D a2 D b2 is SEARCHED so that the same trunk letters c2, c3 progress it too (joint walk, PROOF 3c) and
both splits are clean (lambda(D beta1 h q) <= 4 and a beta2 with D beta2 h q in BK, lambda(D beta2 h p) <= 4).
Then: injections of x1 and x2 into t1 and/or t2 (n_inj random letters each, half sigma3-twisted), n_side side terms at the
ends, and the layered LM d = 5..1 (stop at first failure).
usage: bf7b.py ntrials seed n_inj n_side nq"""
import sys, random, time
from b5eng import *
from lsearch2 import rand_free, dstep
from census import geo
from bf7 import sigma3, mul, solve


def find(E, rng, nq, nlet=4000):
    m = E.m
    b1 = geo(E, rng, E.ID); g1 = E.intern(E.snap(E.D @ b1))
    a1 = geo(E, rng, g1); p = mul(E, a1, g1)
    c2 = geo(E, rng, p); p2 = mul(E, c2, p)
    c3 = geo(E, rng, p2); p3 = mul(E, c3, p2)
    be1 = geo(E, rng, p3)
    if any(v is None for v in (b1, a1, c2, c3, be1)):
        return None, 'no geodesic'
    hits = 0; stat = [0, 0, 0, 0]
    for _ in range(nq):
        b2 = rand_free(E, rng)
        if dstep(E, b2, E.ID) != 1:
            continue
        g2 = E.intern(E.snap(E.D @ b2))
        a2 = rand_free(E, rng)
        if dstep(E, a2, g2) != 1:
            continue
        q = mul(E, a2, g2)
        if q == p:
            continue
        stat[0] += 1
        if dstep(E, c2, q) != 1:
            continue
        stat[1] += 1
        q2 = mul(E, c2, q)
        if dstep(E, c3, q2) != 1:
            continue
        stat[2] += 1
        q3 = mul(E, c3, q2)
        if E.layer(mul(E, be1, q3)) > m - 1:
            continue
        stat[3] += 1
        for _ in range(nlet):
            be2 = rand_free(E, rng)
            if E.dK(mul(E, be2, q3)) == 0:
                if E.layer(mul(E, be2, p3)) <= m - 1:
                    return dict(b1=b1, b2=b2, a1=a1, a2=a2, c2=c2, c3=c3, e1=be1, e2=be2), stat
                break          # clean-ness of the split does not depend on the letter (PROOF 3c): reject this q
    return None, stat


def build(E, rng, sk, n_inj, n_side, where):
    P = Prog(E, 2); sg = sigma3(E)
    s1 = P.add_move([(0, sk['b1'])]); s2 = P.add_move([(1, sk['b2'])])

    def inj():
        out = []
        for l in (0, 1):
            for _ in range(n_inj):
                b = rand_free(E, rng)
                if rng.random() < 0.5:
                    b = sg @ b
                out.append((l, b))
        return out
    t1 = P.add_move([(s1, sk['a1']), (s2, sk['a2'])] + (inj() if 1 in where else []))
    t2 = P.add_move([(t1, sk['c2'])] + (inj() if 2 in where else []))
    t3 = P.add_move([(t2, sk['c3'])])
    for k in (0, 1):
        summ = [(t3, sk['e%d' % (k + 1)])]
        cand = [0, 1, s1, s2, t1, t2]
        for _ in range(n_side):
            summ.append((cand[rng.randrange(len(cand))], rand_free(E, rng)))
        P.add_move(summ)
    return P


def main():
    ntr, seed, n_inj, n_side, nq = [int(a) for a in sys.argv[1:6]]
    rng = random.Random(seed); t0 = time.time(); nfeas = 0; hist = {}
    for tr in range(ntr):
        E = Eng(5, seed=seed * 1000 + tr)
        sk, stat = find(E, rng, nq)
        if sk is None:
            print("trial %d: no special-position partner (q tried with progress starts, c2 ok, c3 ok, clean1) %s (%.0fs)" % (
                tr, stat, time.time() - t0), flush=True)
            continue
        nfeas += 1
        # Lemma Delta check: r = q^-1 p, its layer l predicts failure on layers >= 5 - l for the pure skeleton
        pm = E.MATS[mul(E, sk['a1'], E.intern(E.snap(E.D @ sk['b1'])))]
        qm = E.MATS[mul(E, sk['a2'], E.intern(E.snap(E.D @ sk['b2'])))]
        rr = E.intern(E.snap(np.conj(qm).T @ pm))
        print("trial %d: FOUND top-feasible skeleton, search stats %s; layer(r = q^-1 p) = %d -> Lemma Delta: pure skeleton fails on some layer >= %d" % (
            tr, stat, E.layer(rr), 5 - E.layer(rr)), flush=True)
        for where in ((), (1,), (2,), (1, 2)):
            P = build(E, rng, sk, n_inj if where else 0, n_side, where)
            if not P.setup_final(cap=3000, maxcells=4e7, maxterms=4000, dmin=1):
                print("  inj %s: refused %s" % (where, P.why), flush=True); continue
            nf, res = solve(P, rng)
            deep = min([d for d in res if res[d] < 1e-10], default=6)
            hist[(where, deep)] = hist.get((where, deep), 0) + 1
            print("  inj at %s: unknowns %d terms %s closure %s | %s (%.0fs)" % (
                where, nf, [len(g) for g in P.gens][2:], [(f[0], f[1]) for f in P.fin],
                " ".join("d%d=%.1e" % (d, res[d]) for d in sorted(res, reverse=True)), time.time() - t0), flush=True)
    print("SUMMARY bf7b n_inj %d n_side %d: top-feasible %d of %d; (injection places, deepest layer passed; 6 = none) %s" % (
        n_inj, n_side, nfeas, ntr, sorted(hist.items())), flush=True)


if __name__ == "__main__":
    main()
