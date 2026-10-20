"""bf7.py (mw2-b5beyond, resumed 14 Oct ~15:15): the ONLY butterfly sub-family that the theory leaves for 2 in 7 at block 5,
built in special position and solved layer by layer.

Theory (PROOF.md 3b, 3c):
  - top layer of a single-letter-merge butterfly: the main stream (p, q) must walk JOINTLY (one letter progressing both
    leads) from the merge to distance 1 from both targets, then end with CLEAN split letters (array k reaches K, the
    other coordinate leaves no layer-m term).  Random merges almost never allow two joint steps (census.py: 0 of 10).
  - Lemma Delta with r = q^-1 p and with r' = p^-1 q: a unit needs x2 injected into t1 or t2 AND x1 injected into t1 or t2.
This script: search merges (biased to special position: b2 = b1 and a2 = a1 tau with tau a short free letter, or random)
for a top-feasible skeleton (joint c2, joint c3, clean split e1, e2), then add injections of x1, x2 into t1 and t2
(n_inj random letters each, half sigma-twisted) and n_side side terms at the ends, and run the layered LM d = 5..1.
usage: bf7.py ntrials seed n_inj n_side nlet"""
import sys, random, time
from b5eng import *
from lsearch2 import rand_free, dstep
from census import geo


def sigma3(E):
    Pm = [1 << t for t in range(E.m)]
    Pm[0], Pm[1], Pm[2] = 1 << 1, 1 << 2, 1 << 0
    return Letters(E).addr(Pm)


def mul(E, b, g):
    return E.intern(E.snap(E.D @ b @ E.MATS[g]))


def short_letter(E, rng):
    Lt = Letters(E); m = E.m; k = rng.randrange(4)
    if k == 0:
        return E.Xs[rng.randrange(1, E.n)]                    # Pauli X
    if k == 1:
        return Lt.zp(rng.randrange(1, E.n))                    # Pauli Z
    if k == 2:
        return Lt.sphase(rng.randrange(1, E.n))                # Clifford S-type phase
    lv = levi_menu(m); name = rng.choice(sorted(lv))           # odd-order Levi element
    return Lt.addr(lv[name])


def find_top_skeleton(E, rng, nlet, biased):
    b1 = geo(E, rng, E.ID)
    g1 = E.intern(E.snap(E.D @ b1))
    a1 = geo(E, rng, g1)
    if biased:
        b2 = b1
        for _ in range(200):
            a2 = a1 @ short_letter(E, rng)
            if dstep(E, a2, g1) == 1:
                break
        else:
            return None
        g2 = g1
    else:
        b2 = geo(E, rng, E.ID); g2 = E.intern(E.snap(E.D @ b2)); a2 = geo(E, rng, g2)
    p = mul(E, a1, g1); q = mul(E, a2, g2)
    if p == q:
        return None
    found = None
    seen = set()
    for _ in range(nlet):
        c = rand_free(E, rng)
        if dstep(E, c, p) != 1 or dstep(E, c, q) != 1:
            continue
        p2, q2 = mul(E, c, p), mul(E, c, q)
        if (p2, q2) in seen:
            continue
        seen.add((p2, q2))
        if len(seen) > 25:
            break
        for _ in range(nlet):
            c3 = rand_free(E, rng)
            if dstep(E, c3, p2) != 1 or dstep(E, c3, q2) != 1:
                continue
            p3, q3 = mul(E, c3, p2), mul(E, c3, q2)
            cl = [None, None]
            for _ in range(nlet):
                be = rand_free(E, rng)
                hp, hq = mul(E, be, p3), mul(E, be, q3)
                if cl[0] is None and E.dK(hp) == 0 and E.layer(hq) <= E.m - 1:
                    cl[0] = be
                if cl[1] is None and E.dK(hq) == 0 and E.layer(hp) <= E.m - 1:
                    cl[1] = be
                if cl[0] is not None and cl[1] is not None:
                    return dict(b1=b1, b2=b2, a1=a1, a2=a2, c2=c, c3=c3, e1=cl[0], e2=cl[1])
            break
    return {'joint2': len(seen)}


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


def solve(P, rng):
    fi = P.free_idx(); nf = len(fi); res = {}
    for d in range(5, 0, -1):
        best = np.inf
        for st in range(3):
            x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
            x, cst = P.lm(d, x0, fi, iters=80)
            best = min(best, cst)
            if best < 1e-20:
                break
        res[d] = best
        if best > 1e-10:
            break
    return nf, res


def main():
    ntr, seed, n_inj, n_side, nlet = [int(a) for a in sys.argv[1:6]]
    rng = random.Random(seed); t0 = time.time(); nfeas = 0; hist = {}
    for tr in range(ntr):
        E = Eng(5, seed=seed * 1000 + tr)
        biased = tr % 2 == 0
        sk = find_top_skeleton(E, rng, nlet, biased)
        if sk is None or 'joint2' in sk:
            print("trial %d (%s): no top-feasible skeleton (%s)" % (tr, 'biased' if biased else 'random', sk), flush=True)
            continue
        nfeas += 1
        for where in ((), (1,), (2,), (1, 2)):
            P = build(E, rng, sk, n_inj if where else 0, n_side, where)
            if not P.setup_final(cap=3000, maxcells=4e7, maxterms=4000, dmin=1):
                print("trial %d inj %s: refused %s" % (tr, where, P.why), flush=True); continue
            nf, res = solve(P, rng)
            deep = min([d for d in res if res[d] < 1e-10], default=6)
            hist[(where, deep)] = hist.get((where, deep), 0) + 1
            print("trial %d (%s) inj at %s: unknowns %d terms %s closure %s | %s (%.0fs)" % (
                tr, 'biased' if biased else 'random', where, nf, [len(g) for g in P.gens][2:],
                [(f[0], f[1]) for f in P.fin], " ".join("d%d=%.1e" % (d, res[d]) for d in sorted(res, reverse=True)),
                time.time() - t0), flush=True)
    print("SUMMARY bf7 n_inj %d n_side %d: top-feasible %d of %d; (injection places, deepest layer passed) %s" % (
        n_inj, n_side, nfeas, ntr, sorted(hist.items())), flush=True)


if __name__ == "__main__":
    main()
