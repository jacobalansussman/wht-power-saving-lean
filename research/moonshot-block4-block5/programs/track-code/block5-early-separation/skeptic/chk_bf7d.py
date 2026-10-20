"""chk_bf7d.py (v-mw2-b5beyond, skeptic): re-run mw2-b5beyond's bf7d.py EXACTLY (same seeds, same rng consumption) and
instrument the partner search:
  - print every leftover row (array l, layer, gid) at the failing layer, and whether bf7d SKIPS it because it is
    already in its `tried` set (rows tried in an earlier round are never retried);
  - for every row: the exact existence criterion for an injection partner at (site, end k): layer(W D) == layer(g)
    (B = monomial Cliffords is the stabiliser of L0 and acts transitively on each distance sphere), and for a side term
    from term h: layer(g h^-1) == 1; plus the 20,000-pool result and a BIG independent pool (separate rng, does not
    disturb bf7d's random stream) to confirm existence constructively.
usage: chk_bf7d.py trial_index seed nq npool rounds bigpool"""
import sys, random, time
import numpy as np
from b5eng import *
from lsearch2 import rand_free, dstep
from bf7 import mul, solve
from bf7b import find
from bf7c import word
from bf7d import Prog2, residual_rows, solve_keep


def big_match(E, Wmat, target_mat, n, seed):
    """search n letters from an independent rng (E.rng swapped out and restored)"""
    saved = E.rng; E.rng = random.Random(seed); r2 = random.Random(seed + 1)
    found = None
    try:
        WD = Wmat @ E.D; Ti = np.conj(target_mat).T
        done = 0
        while done < n and found is None:
            nb = min(5000, n - done)
            pool = np.stack([rand_free(E, r2) for _ in range(nb)])
            G = np.einsum('ij,njk->nik', WD, pool) @ Ti
            nzr = (np.abs(G) > 1e-9).sum(axis=2)
            ok = np.nonzero((nzr == 1).all(axis=1))[0]
            for i in ok:
                if E.is_mono(E.snap(G[i])):
                    found = pool[i]; break
            done += nb
    finally:
        E.rng = saved
    return found, done


def main():
    tri, seed, nq, npool, rounds, bigpool = [int(a) for a in sys.argv[1:7]]
    rng = random.Random(seed); t0 = time.time()
    for tr in range(tri + 1):
        E = Eng(5, seed=seed * 1000 + tr)
        sk, stat = find(E, rng, nq)
        if sk is None:
            print("trial %d: no skeleton %s" % (tr, stat), flush=True); continue
        pool = np.stack([rand_free(E, rng) for _ in range(npool)])
        Ws = {(1, 0): word(E, [sk['e1'], sk['c3'], sk['c2']]), (1, 1): word(E, [sk['e2'], sk['c3'], sk['c2']]),
              (2, 0): word(E, [sk['e1'], sk['c3']]), (2, 1): word(E, [sk['e2'], sk['c3']])}
        WDlay = {key: E.layer(E.intern(E.snap(Ws[key] @ E.D))) for key in Ws}

        def match(Wmat, target_mat):
            G = np.einsum('ij,njk->nik', Wmat @ E.D, pool) @ np.conj(target_mat).T
            nzr = (np.abs(G) > 1e-9).sum(axis=2)
            ok = np.nonzero((nzr == 1).all(axis=1))[0]
            for i in ok:
                if E.is_mono(E.snap(G[i])):
                    return pool[i]
            return None
        inj = {1: [], 2: []}; extra = [[], []]; tried = set()
        rowhist = []
        for rd in range(rounds + 1):
            P = Prog2(E, 2)
            s1 = P.add_move([(0, sk['b1'])]); s2 = P.add_move([(1, sk['b2'])])
            t1 = P.add_move([(s1, sk['a1']), (s2, sk['a2'])] + inj[1])
            t2 = P.add_move([(t1, sk['c2'])] + inj[2])
            t3 = P.add_move([(t2, sk['c3'])])
            ends = [P.add_move([(t3, sk['e1'])] + extra[0]), P.add_move([(t3, sk['e2'])] + extra[1])]
            if not P.setup_final():
                print("  round %d: refused %s" % (rd, P.why), flush=True); break
            res, dfail, alpha = solve_keep(P, rng)
            if tr < tri:
                # trials before the one we study: replay exactly as bf7d (same rng use), no instrumentation
                pass
            print("trial %d round %d: inj t1 %d t2 %d, side %s, closure %s | %s (%.0fs)" % (
                tr, rd, len(inj[1]), len(inj[2]), [len(x) for x in extra], [(f[0], f[1]) for f in P.fin],
                " ".join("d%d=%.1e" % (d, res[d]) for d in sorted(res, reverse=True)), time.time() - t0), flush=True)
            if dfail is None or rd == rounds:
                break
            rows = residual_rows(P, alpha, dfail)
            rowhist.append(list(rows))
            if tr == tri:
                print("  WD layers per (site,end): %s" % WDlay, flush=True)
                for (l, g) in rows:
                    lg = E.layer(g); gm = E.MATS[g]
                    crit = {key: (WDlay[key] == lg) for key in Ws}
                    poolhit = {key: (match(Ws[key], gm) is not None) for key in Ws}
                    side = []
                    for i in (s1, s2, t1, t2):
                        for (l2, h) in P.gens[i]:
                            if l2 == l and E.layer(E.intern(E.snap(gm @ np.conj(E.MATS[h]).T))) == 1:
                                side.append(i)
                    print("  ROW l=%d layer=%d gid=%d  already_tried=%s  inj-criterion %s  20k-pool-hit %s  side-term sources %s" % (
                        l, lg, g, (l, g) in tried, crit, poolhit, side), flush=True)
                    for key in Ws:
                        if crit[key] and not poolhit[key]:
                            b, n = big_match(E, Ws[key], gm, bigpool, 7 + g)
                            print("     big pool for %s: %s after %d letters" % (key, 'FOUND' if b is not None else 'none', n), flush=True)
            # ---- bf7d's partner search, verbatim (consumes no rng) ----
            added = 0
            for (l, g) in rows[:12]:
                if (l, g) in tried:
                    continue
                tried.add((l, g)); gm = E.MATS[g]
                done = False
                for (site, k) in ((2, 0), (2, 1), (1, 0), (1, 1)):
                    b = match(Ws[(site, k)], gm)
                    if b is not None:
                        inj[site].append((l, b)); added += 1; done = True; break
                if done:
                    continue
                for k in (0, 1):
                    for i in (s1, s2, t1, t2):
                        for (l2, h) in P.gens[i]:
                            if l2 != l:
                                continue
                            b = match(E.I, gm @ np.conj(E.MATS[h]).T)
                            if b is not None:
                                extra[k].append((i, b)); added += 1; done = True; break
                        if done:
                            break
                    if done:
                        break
            print("  round %d: %d leftover rows at layer %d, %d partners added" % (rd, len(rows), dfail, added), flush=True)
            if tr == tri and len(rowhist) >= 2:
                print("  rows this round identical to previous round's rows: %s" % (set(rowhist[-1]) == set(rowhist[-2])), flush=True)
            if added == 0:
                break
        if tr == tri:
            break
    print("DONE (%.0fs)" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
