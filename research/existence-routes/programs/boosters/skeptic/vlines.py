# v-e3-boost vlines.py: NEW certificate family, permutation module on ALL lines (2-dim subspaces) of the label space
# F_2^(2m), Y_U = span{e_W : type_U(W) in T}, type_U(W) = (dim W cap U, omega|_W) in {(0,0),(0,1),(1,0),(1,1),(2,0)}.
# Valid by Theorem R: free letters fix L_0 and are symplectic, so they keep Y_{L_0}; f(g) = |O(gL_0) \ O(L_0)|, which
# depends only on d(L_0, gL_0) (Sp(2m,2) is transitive on frames at each distance from L_0; O is Sp-equivariant).
# Counting by a row-by-row DP over ordered pairs (a, b) relative to U = L_0 (Z-type) and V_d = H on rows 1..d of L_0.
import sys
from itertools import product, combinations
TYPES = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0)]
def counts(m, d):
    # per row i: a_i = (ax, az), b_i = (bx, bz) in F_2^2 each.  U: x = 0.  V_d: rows <= d: z = 0; rows > d: x = 0.
    # state: (a nz, b nz, a!=b, for c in (a, b, a+b): inU, inV; omega parity)
    from collections import defaultdict
    st = defaultdict(int); st[(0, 0, 0, (1, 1, 1), (1, 1, 1), 0)] = 1
    for i in range(m):
        new = defaultdict(int)
        for (anz, bnz, neq, inU, inV, om), cnt in st.items():
            for ax, az, bx, bz in product((0, 1), repeat=4):
                cx, cz = ax ^ bx, az ^ bz
                vecs = ((ax, az), (bx, bz), (cx, cz))
                nU = tuple(inU[t] & (vecs[t][0] == 0) for t in range(3))
                if i < d: nV = tuple(inV[t] & (vecs[t][1] == 0) for t in range(3))
                else:     nV = tuple(inV[t] & (vecs[t][0] == 0) for t in range(3))
                key = (anz | ax | az, bnz | bx | bz, neq | cx | cz, nU, nV, om ^ ((ax & bz) ^ (az & bx)))
                new[key] += cnt
        st = new
    J = {}
    for (anz, bnz, neq, inU, inV, om), cnt in st.items():
        if not (anz and bnz and neq): continue
        tU = (sum(inU) // 1, om); tV = (sum(inV), om)
        jU = {0: 0, 1: 1, 3: 2}[sum(inU)]; jV = {0: 0, 1: 1, 3: 2}[sum(inV)]
        key = ((jU, om), (jV, om))
        J[key] = J.get(key, 0) + cnt
    return {k: v // 6 for k, v in J.items()}   # each line = 6 ordered pairs
def profile(m, T):
    out = []
    for d in range(1, m + 1):
        J = counts(m, d)
        inV = sum(c for (tu, tv), c in J.items() if tv in T)
        both = sum(c for (tu, tv), c in J.items() if tv in T and tu in T)
        out.append(inV - both)
    return out
if __name__ == "__main__":
    # sanity: total lines = (4^m-1)(4^m-2)/6; d = 0 (V = U) gives zero profile
    for m in (2, 3):
        J = counts(m, 0); tot = sum(J.values())
        print(f"m={m}: lines {tot} expect {(4**m-1)*(4**m-2)//6}; diagonal only at d=0: {all(a == b for (a, b) in J)}")
    # brute-force check at m = 3 for every T
    def brute(m, d, T):
        n = 2 * m
        def inU(v): return (v >> m) == 0
        mask = (1 << d) - 1
        def inV(v): x, z = v >> m, v & ((1 << m) - 1); return (z & mask) == 0 and (x & ~mask) == 0
        def om(a, b):
            ax, az, bx, bz = a >> m, a & ((1 << m) - 1), b >> m, b & ((1 << m) - 1)
            return (bin(ax & bz).count("1") + bin(az & bx).count("1")) & 1
        lines = set()
        for a in range(1, 1 << n):
            for b in range(a + 1, 1 << n):
                lines.add(frozenset((a, b, a ^ b)))
        cU = cB = 0
        for L in lines:
            a, b, _ = sorted(L)
            tU = ({0: 0, 1: 1, 3: 2}[sum(map(inU, L))], om(a, b)); tV = ({0: 0, 1: 1, 3: 2}[sum(map(inV, L))], om(a, b))
            if tV in T:
                cU += 1
                if tU in T: cB += 1
        return cU - cB
    ok = True
    for r in range(1, 6):
        for T in combinations(TYPES, r):
            for d in (1, 2, 3):
                if brute(3, d, set(T)) != profile(3, set(T))[d - 1]: ok = False
    print("m=3 brute force = DP for every T and d:", ok)
