"""bf7r.py (mw2-b5beyond, resumed 14 Oct ~15:15): COLLISION REPAIR on special-position butterflies (bf7b.py skeletons).
After the layered solve fails at layer d, every term g0 of an end e_k with layer >= d (other than its K-aligned term) is
offered a partner: an earlier generator term h (same coordinate) with D^-1 g0 h^-1 MONOMIAL (a free letter b); the side
term (that generator, b) is added to e_k's content so that the coefficients can cancel g0 (a collision, rule C1).  Up to
`rounds` rounds (new side terms bring new terms).  Injections of x1, x2 into t1/t2 as in bf7b (B++ requires them).
Lemma Delta predicts: with NO injection the repair cannot pass layer 3 (end-only repairs carry defect <= 4).
usage: bf7r.py ntrials seed n_inj nq rounds"""
import sys, random, time
from b5eng import *
from lsearch2 import rand_free, dstep
from census import geo
from bf7 import sigma3, mul, solve
from bf7b import find


def build(E, rng, sk, n_inj, where, extra):
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
    t1 = P.add_move([(s1, sk['a1']), (s2, sk['a2'])] + (sk.setdefault('inj1', inj()) if 1 in where else []))
    t2 = P.add_move([(t1, sk['c2'])] + (sk.setdefault('inj2', inj()) if 2 in where else []))
    t3 = P.add_move([(t2, sk['c3'])])
    ends = []
    for k in (0, 1):
        ends.append(P.add_move([(t3, sk['e%d' % (k + 1)])] + extra[k]))
    return P, ends


def partners(E, P, ends, d):
    """side terms (src, letter) for each end that collide with the end's own terms of layer >= d (not K-aligned)."""
    D_inv = np.conj(E.D).T
    out = [[], []]
    for k, e in enumerate(ends):
        for (l, g0) in P.gens[e]:
            if E.layer(g0) < d or E.dK(g0) == 0:
                continue
            for i in range(e):
                if i in ends:
                    continue
                for (l2, h) in P.gens[i]:
                    if l2 != l:
                        continue
                    M = E.snap(D_inv @ E.MATS[g0] @ np.conj(E.MATS[h]).T)
                    if E.is_mono(M):
                        out[k].append((i, M))
    return out


def main():
    ntr, seed, n_inj, nq, rounds = [int(a) for a in sys.argv[1:6]]
    rng = random.Random(seed); t0 = time.time(); hist = {}
    for tr in range(ntr):
        E = Eng(5, seed=seed * 1000 + tr)
        sk, stat = find(E, rng, nq)
        if sk is None:
            print("trial %d: no special-position partner %s (%.0fs)" % (tr, stat, time.time() - t0), flush=True)
            continue
        print("trial %d: skeleton found %s" % (tr, stat), flush=True)
        for where in ((), (1, 2)):
            extra = [[], []]; seen = [set(), set()]; last = None
            for rd in range(rounds + 1):
                P, ends = build(E, rng, sk, n_inj, where, extra)
                if not P.setup_final(cap=3000, maxcells=4e7, maxterms=6000, dmin=1):
                    print("  inj %s round %d: refused %s" % (where, rd, P.why), flush=True); break
                nf, res = solve(P, rng)
                deep = min([d for d in res if res[d] < 1e-10], default=6)
                fail = min(res)
                print("  inj %s round %d: side terms %s unknowns %d closure %s | %s (%.0fs)" % (
                    where, rd, [len(x) for x in extra], nf, [(f[0], f[1]) for f in P.fin],
                    " ".join("d%d=%.1e" % (d, res[d]) for d in sorted(res, reverse=True)), time.time() - t0), flush=True)
                last = deep
                if deep == 1 or rd == rounds:
                    break
                new = partners(E, P, ends, fail)
                added = 0
                for k in (0, 1):
                    for (i, M) in new[k]:
                        key = (i, E.key(M))
                        if key in seen[k] or len(extra[k]) >= 12:
                            continue
                        seen[k].add(key); extra[k].append((i, M)); added += 1
                if added == 0:
                    print("  inj %s round %d: no colliding partner for the leftovers at layer %d" % (where, rd, fail), flush=True)
                    break
            hist[(where, last)] = hist.get((where, last), 0) + 1
    print("SUMMARY bf7r: (injection places, deepest layer passed after repair; 6 = none) %s" % sorted(hist.items()), flush=True)


if __name__ == "__main__":
    main()
