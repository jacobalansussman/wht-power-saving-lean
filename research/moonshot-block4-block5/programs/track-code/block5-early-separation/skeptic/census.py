"""census.py (mw2-b5beyond, resumed 14 Oct ~15:05): frame geometry of a SHARED single-stream trunk.
After a single-letter merge t1 = D(a1 s1 + a2 s2) the register is (p, q) and every later pure-trunk register is (y p, y q):
the two arrays move rigidly (relative factor r = q^-1 p fixed).  The top layer then needs, at the ends, CLEAN SPLIT letters
beta: D beta g_k in B K (array k reaches K) while lambda(D beta g_k') <= m-1 (array k' leaves NO top-layer term), because a
top term of the other coordinate is rigidly attached and cannot be cancelled (see PROOF.md section 3c).
For random geodesic starts and merges (a_j progress for s_j's lead) this script samples N random free letters per step and counts:
  joint  : letters progressing both leads (dK - 1 for both)          -- needed at t2, t3 (single-letter trunk)
  clean_k: at the trunk end, letters with D beta g_k in BK and lambda(D beta g_k') <= m-1
over block m and trunk length (m-3 joint steps between the merge and the ends: 2 in m+2 butterfly).
usage: census.py m ntrials nletters seed"""
import sys, random
from b5eng import *
from lsearch2 import rand_free, dstep


def geo(E, rng, g, tries=4000):
    for _ in range(tries):
        b = rand_free(E, rng)
        if dstep(E, b, g) == 1:
            return b
    return None


def main():
    m, ntr, nlet, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    rng = random.Random(seed)
    stats = {'trials': 0, 'joint_ok': 0, 'clean_both': 0, 'clean_one': 0}
    for tr in range(ntr):
        E = Eng(m, seed=seed * 1000 + tr)
        ID = E.ID
        b1 = geo(E, rng, ID); b2 = geo(E, rng, ID)
        g1 = E.intern(E.snap(E.D @ b1)); g2 = E.intern(E.snap(E.D @ b2))
        a1 = geo(E, rng, g1); a2 = geo(E, rng, g2)
        p = E.intern(E.snap(E.D @ a1 @ E.MATS[g1])); q = E.intern(E.snap(E.D @ a2 @ E.MATS[g2]))
        lay = (E.layer(p), E.layer(q)); dk = (E.dK(p), E.dK(q))
        line = ["trial %d m=%d merge leads layer %s dK %s" % (tr, m, lay, dk)]
        stats['trials'] += 1
        ok = True
        for step in range(m - 3):        # joint single-letter steps t2 .. t_{m-2}
            nj = 0; first = None
            for _ in range(nlet):
                b = rand_free(E, rng)
                if dstep(E, b, p) == 1 and dstep(E, b, q) == 1:
                    nj += 1
                    if first is None:
                        first = b
            line.append("step %d joint %d/%d" % (step + 2, nj, nlet))
            if first is None:
                ok = False; break
            p = E.intern(E.snap(E.D @ first @ E.MATS[p])); q = E.intern(E.snap(E.D @ first @ E.MATS[q]))
        if ok:
            stats['joint_ok'] += 1
            cl = [0, 0]; geo_k = [0, 0]
            for _ in range(nlet):
                b = rand_free(E, rng)
                hp = E.intern(E.snap(E.D @ b @ E.MATS[p])); hq = E.intern(E.snap(E.D @ b @ E.MATS[q]))
                if E.dK(hp) == 0:
                    geo_k[0] += 1
                    if E.layer(hq) <= m - 1:
                        cl[0] += 1
                if E.dK(hq) == 0:
                    geo_k[1] += 1
                    if E.layer(hp) <= m - 1:
                        cl[1] += 1
            line.append("trunk-end leads layer %s dK %s; reach-K letters %s, CLEAN split %s of %d" % (
                (E.layer(p), E.layer(q)), (E.dK(p), E.dK(q)), geo_k, cl, nlet))
            if cl[0] and cl[1]:
                stats['clean_both'] += 1
            elif cl[0] or cl[1]:
                stats['clean_one'] += 1
        print(" | ".join(line), flush=True)
    print("SUMMARY m=%d: %s" % (m, stats), flush=True)


if __name__ == "__main__":
    main()
