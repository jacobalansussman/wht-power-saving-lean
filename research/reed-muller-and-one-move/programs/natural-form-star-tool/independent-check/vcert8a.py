"""vcert8a.py: line-by-line check of a starrec2 CERT clique/matching certificate with own code (vf + vrec flats).
For each cited cell x,y:F/S/round/d/i : r_F = nu (<= 0 here) for the target set; S = D_F n A for a d-flat A;
own round-1 recursion of F alone gives u(A) = 0 (or, d = -1, the quadratic star with G1, G2 gives 0); (x,y) in SxS.
Then the clique: all |Y|^2 cells cited, |Y| > W."""
import sys, numpy as np, vf, vrec
def rec_u(F, w, nu):
    P = vrec.pre(F); u = []
    for dd in range(5):
        mw, sz, S = P[dd]; u.append(np.minimum(sz, (mw - w[None, :]).min(axis=1)).astype(np.int64))
    for dd in range(3, -1, -1):
        g = (1 << (5 - dd)) - 1; sz = P[dd][1]; Bu = u[dd + 1][vrec.UP[dd]]
        val = Bu.sum(axis=1) - 2 * sz + 2 * nu
        assert not ((val < 0) & (sz > 0)).any()
        u[dd] = np.where(sz > 0, np.minimum(u[dd], np.minimum(val // (g - 2), Bu.min(axis=1))), u[dd])
    return u, P
def main(line):
    t = line.split(); ks = [int(x) for x in t[1].split(',')]; assert t[2] == 'C'
    W = len(ks); q = int(t[4]); Y = int(t[6][2:]); Yd = [z for z in range(1, 32) if Y >> (z - 1) & 1]
    w = vrec.CLSPAR[:, ks].sum(axis=1); r = w[vrec.VID]
    assert len(Yd) == q and q > W
    cells = {}; cache = {}; okc = 0
    for c in t[7:]:
        xy, rest = c.split(':'); x, y = map(int, xy.split(',')); F, S, rnd, d, i = map(int, rest.split('/'))
        assert rnd == 1
        key = (F, S, d)
        if key not in cache:
            nu = int(r[F]); assert nu == 0, ("test form not even", F)
            assert S & ~int(vf.MASK[F]) == 0
            if d >= 0:
                u, P = cache.get(F) or rec_u(F, w, nu)
                hits = [j for j in range(len(P[d][2])) if int(P[d][2][j]) == S and u[d][j] == 0]
                cache[key] = len(hits) > 0
            else:
                cache[key] = None
        assert cache[key], ("zero block not re-derived", F, S, d)
        assert (S >> (x - 1) & 1) and (S >> (y - 1) & 1)
        cells[(x, y)] = 1; okc += 1
    allc = all((x, y) in cells for x in Yd for y in Yd)
    print("set", t[1], "clique", q, ">", W, "| cited cells", okc, "| all |Y|^2 cells cited:", allc, "| distinct blocks re-derived:", sum(1 for v in cache.values() if v))
    return allc
if __name__ == "__main__":
    for l in open(sys.argv[1]):
        if l.startswith('S ') and ' C ' in l: main(l)
