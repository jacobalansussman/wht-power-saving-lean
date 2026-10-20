"""bf7c.py (mw2-b5beyond, resumed 14 Oct ~15:15): special-position butterflies (bf7b) with COLLISION-DESIGNED injections.
The clean split leaves, in end e_1, the coordinate-2 term L_2 = D beta_1 h q (layer 3 or 4), rigidly tied to array 1's K term,
and in e_2 the coordinate-1 term L_1 = D beta_2 h p.  Output 1 must cancel L_2 with a term of ANOTHER stream in the same
left coset B L_2 (PROOF 3c').  An x_2 stream injected at t_1 with letter beta reaches e_1 as W D beta (W = D beta_1 D c_3 D c_2):
it collides with L_2 iff frame(W D beta) = frame(L_2), possible only if d(D L_0, F_W) = d(q L_0, F_W) (same B-orbit).
This script searches beta (and the analogous x_1 letter for L_1, and injections at t_2: W = D beta_1 D c_3) among random
free letters by FRAME MATCHING, adds the colliding injections, and runs the layered LM d = 5..1.
usage: bf7c.py ntrials seed nq nbeta"""
import sys, random, time
from b5eng import *
from lsearch2 import rand_free, dstep
from census import geo
from bf7 import mul, solve
from bf7b import find


def word(E, letters):
    """D l_1 D l_2 ... as a matrix (applied right to left: last letter first)"""
    M = E.I
    for b in letters:
        M = E.snap(M @ E.D @ b)
    return M


def match(E, rng, W, target, nbeta):
    """letters b with W D b in the left coset B target  <=>  (W D b) target^-1 monomial (Z4 entries)."""
    Ti = np.conj(E.MATS[target]).T; WD = W @ E.D; hits = []
    for _ in range(nbeta):
        b = rand_free(E, rng)
        G = WD @ b @ Ti
        nz = np.abs(G) > 1e-9
        if (nz.sum(axis=1) == 1).all() and E.is_mono(E.snap(G)):
            hits.append(b)
            if len(hits) >= 2:
                break
    return hits


def main():
    ntr, seed, nq, nbeta = [int(a) for a in sys.argv[1:5]]
    rng = random.Random(seed); t0 = time.time(); hist = {}
    for tr in range(ntr):
        E = Eng(5, seed=seed * 1000 + tr)
        sk, stat = find(E, rng, nq)
        if sk is None:
            print("trial %d: no skeleton %s" % (tr, stat), flush=True); continue
        g1 = E.intern(E.snap(E.D @ sk['b1'])); g2 = E.intern(E.snap(E.D @ sk['b2']))
        p = mul(E, sk['a1'], g1); q = mul(E, sk['a2'], g2)
        h = word(E, [sk['c3'], sk['c2']])
        L2 = E.intern(E.snap(E.D @ sk['e1'] @ h @ E.MATS[q]))
        L1 = E.intern(E.snap(E.D @ sk['e2'] @ h @ E.MATS[p]))
        W1a = word(E, [sk['e1'], sk['c3'], sk['c2']]); W2a = word(E, [sk['e2'], sk['c3'], sk['c2']])
        W1b = word(E, [sk['e1'], sk['c3']]); W2b = word(E, [sk['e2'], sk['c3']])
        inj = {1: [], 2: []}; found = []
        for (site, W, L, lin, name) in ((1, W1a, L2, 1, 'x2@t1->L2'), (1, W2a, L1, 0, 'x1@t1->L1'),
                                        (2, W1b, L2, 1, 'x2@t2->L2'), (2, W2b, L1, 0, 'x1@t2->L1')):
            hs = match(E, rng, W, L, nbeta)
            if hs:
                inj[site].append((lin, hs[0])); found.append(name)
        print("trial %d: skeleton %s; leftover layers L2 %d, L1 %d; colliding injections found: %s (%.0fs)" % (
            tr, stat, E.layer(L2), E.layer(L1), found, time.time() - t0), flush=True)
        if not found:
            continue
        P = Prog(E, 2)
        s1 = P.add_move([(0, sk['b1'])]); s2 = P.add_move([(1, sk['b2'])])
        t1 = P.add_move([(s1, sk['a1']), (s2, sk['a2'])] + inj[1])
        t2 = P.add_move([(t1, sk['c2'])] + inj[2])
        t3 = P.add_move([(t2, sk['c3'])])
        for k in (0, 1):
            P.add_move([(t3, sk['e%d' % (k + 1)])])
        if not P.setup_final(cap=3000, maxcells=4e7, maxterms=6000, dmin=1):
            print("  refused %s" % P.why, flush=True); continue
        nf, res = solve(P, rng)
        deep = min([d for d in res if res[d] < 1e-10], default=6)
        hist[deep] = hist.get(deep, 0) + 1
        print("  with %s: unknowns %d terms %s closure %s | %s (%.0fs)" % (
            found, nf, [len(g) for g in P.gens][2:], [(f[0], f[1]) for f in P.fin],
            " ".join("d%d=%.1e" % (d, res[d]) for d in sorted(res, reverse=True)), time.time() - t0), flush=True)
    print("SUMMARY bf7c: deepest layer passed with collision-designed injections (6 = none): %s" % sorted(hist.items()), flush=True)


if __name__ == "__main__":
    main()
