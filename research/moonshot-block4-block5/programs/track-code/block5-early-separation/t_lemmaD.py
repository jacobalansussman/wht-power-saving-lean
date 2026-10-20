"""t_lemmaD.py (mw2-b5beyond): numerical check of Corollary B+ (PROOF.md section 3) at block 5.
Butterfly with single-letter starts and merge, a PURE trunk with multi-term letters, ends with many side terms, inputs and
multi-term letters.  Lemma Delta (r of layer <= 2) says the identity restricted to layers >= 3 can never hold
(res_3 > 0 for every coefficient choice); the top layer alone may pass.  Positive control: the same skeleton with ONE input injected into the trunk (t1 gets gamma x_2): the lemma no
longer applies (no claim).  Letters: random free (affine x Clifford phase x Pauli) and geodesic letters for array 1/2 so
that K's frame IS reached (otherwise the test is vacuous).
usage: t_lemmaD.py nskel seed [inject]"""
import sys, random, time
from b5eng import *
from lsearch2 import rand_free, dstep


def geo(E, rng, g, tries=3000):
    for _ in range(tries):
        b = rand_free(E, rng)
        if dstep(E, b, g) == 1:
            return b
    return None


def main():
    nskel, seed = int(sys.argv[1]), int(sys.argv[2]); inject = len(sys.argv) > 3 and sys.argv[3] == '1'
    rng = random.Random(seed); t0 = time.time()
    for s in range(nskel):
        E = Eng(5, seed=seed * 100 + s)
        P = Prog(E, 2)
        # starts: geodesic first letters
        b1 = geo(E, rng, E.ID); b2 = geo(E, rng, E.ID)
        s1 = P.add_move([(0, b1)]); s2 = P.add_move([(1, b2)])
        g1 = P.gens[s1][0][1]; g2 = P.gens[s2][0][1]
        a1 = geo(E, rng, g1); a2 = geo(E, rng, g2)
        summ = [(s1, a1), (s2, a2)]
        if inject:
            summ.append((1, rand_free(E, rng)))
        t1 = P.add_move(summ)
        prev = t1; ok = True
        for k in range(2):           # t2, t3: multi-term letters: one geodesic letter per array lead + one random
            L = {}
            for (l, g) in P.gens[prev]:
                if l not in L or E.dK(g) < E.dK(L[l]):
                    L[l] = g
            bs = [geo(E, rng, L[0]), geo(E, rng, L[1]), rand_free(E, rng)]
            if any(b is None for b in bs):
                ok = False; break
            prev = P.add_move([(prev, b) for b in bs])
        if not ok:
            print("skel %d: no geodesic letter" % s); continue
        t3 = prev
        ends = []
        for j in range(2):
            L = {}
            for (l, g) in P.gens[t3]:
                if l not in L or E.dK(g) < E.dK(L[l]):
                    L[l] = g
            summ = [(t3, geo(E, rng, L[0])), (t3, geo(E, rng, L[1]))]
            for src in [s1, s2, t1, t3 - 1, 0, 1]:
                summ.append((src, rand_free(E, rng)))
            ends.append(P.add_move(summ))
        if not P.setup_final(cap=3000, dmin=3):
            print("skel %d: refused %s" % (s, P.why)); continue
        fi = P.free_idx(); nf = len(fi)
        out = []
        for d in (5, 4, 3):
            best = np.inf
            for st in range(4):
                x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
                x, c = P.lm(d, x0, fi, iters=80)
                best = min(best, c)
            out.append("res_%d=%.2e" % (d, best))
            if best > 1e-10:
                break
        print("skel %d: inject=%d unknowns %d terms %s closure(>=3) %s | %s  (%.0fs)" % (
            s, inject, nf, [len(g) for g in P.gens][2:], [(f[0], f[1]) for f in P.fin], " ".join(out), time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
