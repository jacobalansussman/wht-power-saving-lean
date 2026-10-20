# dpg.py -- the LP on the WHOLE dual polar graph at block 3 (pf-duality, 14 Oct 2026).
# Frames: all 135 Lagrangians of F_2^6 (symplectic form x.z' + x'.z); d(L,L') = 3 - dim(L & L').
# For a certificate (rho, N), Y_L = rho(gamma)N for any Clifford gamma with gamma L_0 = L is well defined (Stab(L_0) keeps N),
# the family is Sp(6,2)-equivariant and r(S) = dim sum Y_L is Sp(6,2)-invariant; delta depends only on the distance
# (PROOF.md Theorem I).  So the unknown r lives on Sp(6,2)-orbits of sets of <= K frames.  Rows as in lpm.py.
import sys, os, itertools, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from math import comb
from lpm import Interval, DFZ, ING, inst_row, solve, FANO_LAB

def omega(u, v):
    return (bin((u & 7) & (v >> 3)).count('1') + bin((v & 7) & (u >> 3)).count('1')) & 1

class DPG(Interval):
    def __init__(self, K, log=print):
        t0 = time.time()
        self.m, self.K = 3, K
        subs = set()
        for a, b, c in itertools.combinations(range(1, 64), 3):
            S = frozenset({0, a, b, c, a ^ b, a ^ c, b ^ c, a ^ b ^ c})
            if len(S) == 8 and all(omega(u, v) == 0 for u in S for v in S): subs.add(S)
        SUBS = sorted(subs, key=lambda S: sorted(S)); NF = len(SUBS); assert NF == 135
        IDX = {S: i for i, S in enumerate(SUBS)}
        self.SUBS, self.NF, self.IDX = SUBS, NF, IDX
        self.dim = np.full(NF, 3)
        D = np.zeros((NF, NF), dtype=np.int64)
        for i in range(NF):
            for j in range(NF): D[i, j] = 3 - (len(SUBS[i] & SUBS[j]).bit_length() - 1)
        self.D = D
        # generators: symplectic transvections u -> u + omega(u,v) v, as permutations of the 135 frames
        def trans(v):
            return [IDX[frozenset(u ^ (v if omega(u, v) else 0) for u in S)] for S in SUBS]
        gens = [np.array(trans(v), dtype=np.int16) for v in (1, 8, 10, 3, 24, 6, 48, 12)]
        ident = np.arange(NF, dtype=np.int16)
        seen = {ident.tobytes()}; frontier = [ident]; allp = [ident]
        while frontier:
            nxt = []
            for p in frontier:
                for g in gens:
                    q = g[p]
                    kb = q.tobytes()
                    if kb not in seen: seen.add(kb); nxt.append(q); allp.append(q)
            frontier = nxt
        self.P = np.array(allp, dtype=np.int16); del seen
        log(f"DPG(3) frames={NF} group={len(self.P)} ({time.time()-t0:.1f}s)")
        assert len(self.P) == 1451520, len(self.P)
        self.Cb = np.array([[comb(a, j) for j in range(K + 2)] for a in range(NF + 1)], dtype=np.int64)
        self.oid, self.reps, self.off = [], [], []
        nv = 0
        for k in range(K + 1):
            N = comb(NF, k); oid = np.full(N, -1, dtype=np.int32); reps = []
            self.off.append(nv)
            if k == 0:
                oid[0] = 0; reps.append(()); nv += 1
            else:
                ptr = 0; CH = 1 << 22
                while True:
                    while ptr < N:
                        nz = np.flatnonzero(oid[ptr:ptr + CH] < 0)
                        if len(nz): ptr += int(nz[0]); break
                        ptr += CH
                    if ptr >= N: break
                    S = self.unrank(ptr, k)
                    imgs = np.sort(self.P[:, S].astype(np.int64), axis=1)
                    oid[self.rank(imgs)] = nv - self.off[k]
                    reps.append(tuple(S)); nv += 1
            self.oid.append(oid); self.reps.append(reps)
            log(f"  size {k}: {len(reps)} orbits ({time.time()-t0:.1f}s)")
        self.NV = nv
        self.EMPTY = 0
    def frame_of(self, vecs):
        S = {0}
        for v in vecs: S |= {x ^ v for x in S}
        return self.IDX[frozenset(S)]

if __name__ == "__main__":
    from cutsearch import run
    K, uniform, rounds, starts, seed, tmax = int(sys.argv[1]), sys.argv[2] == 'uniform', int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    I = DPG(K)
    L0 = I.frame_of([8, 16, 32])          # Z-type: {(0|z)}
    def N(v):                              # interval frame N_<v> = span{(v|0)} + {(0|z): z.v = 0}
        zs = [z for z in range(1, 8) if bin(z & v).count('1') % 2 == 0]
        return I.frame_of([v] + [z << 3 for z in zs])
    letters = [(N(lab),) for lab in FANO_LAB]
    cuts = []
    for nmT, T in (("none", ()), ("L0", (L0,))):
        r = inst_row(I, DFZ, letters, T)
        if r: cuts.append((f"dfz:fano-neighbours@L0,T={nmT}", r))
    out, cuts, hist = run(I, cuts, uniform, rounds, starts, seed, tmax, tag=f"dpg3K{K}")
    json.dump({"K": K, "uniform": uniform, "history": hist, "opt": out["opt"], "delta": out["delta"], "active": out["active"],
               "cuts": [(nm, {str(k): v for k, v in r.items()}) for nm, r in cuts]},
              open(f"cuts_dpg3_K{K}_{'u' if uniform else 'a'}_{seed}.json", "w"))
