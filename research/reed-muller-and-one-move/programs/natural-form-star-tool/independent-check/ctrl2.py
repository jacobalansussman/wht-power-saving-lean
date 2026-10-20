"""ctrl2.py: all frames (5 directions z with sum z z^T = T_k) by brute force; disjoint frame packings for the plain
controls; TRUE E of each plain program; T10 on all forms; vrec.run with trueE (u(A) >= true rank, zeros true)."""
import sys, itertools, numpy as np
import vf, vrec
PAIRS = vf.PAIRS
def tvec_z(z): return sum(((z >> j) & (z >> k) & 1) << p for p, (j, k) in enumerate(PAIRS))
def tvec_T(T): return sum(T[j][k] << p for p, (j, k) in enumerate(PAIRS))
TV = {z: tvec_z(z) for z in range(1, 32)}
def frames(k):
    t = tvec_T(vf.Tmat(k)); out = []
    for c in itertools.combinations(range(1, 32), 5):
        x = 0
        for z in c: x ^= TV[z]
        if x == t: out.append(c)
    return out
def pack(ks):
    F = [frames(k) for k in ks]
    def rec(i, used, acc):
        if i == len(ks): return acc
        for f in F[i]:
            if not used & set(f):
                r = rec(i + 1, used | set(f), acc + [f])
                if r: return r
        return None
    return rec(0, set(), []), [len(f) for f in F]
def trueE(frs, seed=0):
    rng = np.random.default_rng(seed); E = np.eye(31)
    for fr in frs:
        fr = list(fr); rng.shuffle(fr)
        for a, x in enumerate(fr):
            E[x - 1, x - 1] = 0
            for b, y in enumerate(fr):
                if b < a: E[x - 1, y - 1] = 1
                if b > a: E[x - 1, y - 1] = -1
    return E
def T10(ks, E):
    w = vf.PAR[:, ks].sum(axis=1); bad = 0
    for Q in range(1, vf.NF):
        D = [z - 1 for z in range(1, 32) if int(vf.MASK[Q]) >> (z - 1) & 1]
        if np.linalg.matrix_rank(E[np.ix_(D, D)]) != len(D) - w[Q]: bad += 1
    return bad
if __name__ == "__main__":
    maxd = int(sys.argv[1])
    for ks in [[2, 3, 20, 30], [4, 11, 15, 24, 29], [0, 3, 6]]:
        frs, nfr = pack(ks); print(ks, "frames per target", nfr, "packing", frs, flush=True)
        if not frs: continue
        E = trueE(frs)
        print(ks, "T10 mismatches", T10(ks, E), flush=True)
        print(ks, "vrec(true E)", vrec.run(ks, maxd=maxd, rounds=3, trueE=E), flush=True)
        print(ks, "vrec(true E, all r=0 forms)", vrec.run(ks, maxd=maxd, rounds=3, allzero=True, trueE=E), flush=True)
