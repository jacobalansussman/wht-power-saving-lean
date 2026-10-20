"""congr_classes.py (rm2-phase): check by orbit enumeration that the nonsingular symmetric m x m matrices over F_2
fall into congruence classes S -> A S A^T as follows: one class (that of I) for odd m; two classes (I and the
alternating ones) for even m. Generators of GL(m,2): the elementary transvections E_ij (row i += row j)."""
import itertools
def run(m):
    pairs = [(i, j) for i in range(m) for j in range(i, m)]
    def tomat(bits):
        S = [[0]*m for _ in range(m)]
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1: S[i][j] = S[j][i] = 1
        return S
    def key(S): return tuple(tuple(r) for r in S)
    def rank(S):
        rows = [sum(S[i][j] << j for j in range(m)) for i in range(m)]; b = []
        for v in rows:
            for x in b: v = min(v, v ^ x)
            if v: b.append(v)
        return len(b)
    allns = set()
    for bits in range(1 << len(pairs)):
        S = tomat(bits)
        if rank(S) == m: allns.add(key(S))
    def act(S, i, j):      # E S E^T with E = I + e_i e_j^T: row i += row j, then column i += column j
        T = [r[:] for r in S]
        T[i] = [T[i][k] ^ T[j][k] for k in range(m)]
        for r in range(m): T[r][i] ^= T[r][j]
        return T
    seen = set(); classes = []
    for s0 in sorted(allns):
        if s0 in seen: continue
        orb = {s0}; stack = [s0]
        while stack:
            S = [list(r) for r in stack.pop()]
            for i in range(m):
                for j in range(m):
                    if i != j:
                        k = key(act(S, i, j))
                        if k not in orb: orb.add(k); stack.append(k)
        seen |= orb
        alt = all(s0[i][i] == 0 for i in range(m))
        classes.append((len(orb), "alternating" if alt else "non-alternating"))
    print(f"m={m}: {len(allns)} nonsingular symmetric; congruence classes {classes}")
for m in (2, 3, 4, 5):
    run(m)
