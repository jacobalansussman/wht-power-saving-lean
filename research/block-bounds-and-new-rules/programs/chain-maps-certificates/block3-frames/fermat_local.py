"""xc3-cert3: how many of its 15 Fermat neighbours can an 'other' frame U be clean with, if the 405 Fermat lifts are
kept?  (mod p = 2^31 - 1, exploration.)

If U is clean with Fermat neighbours V_j (j in J) then W_j = Y_U cap Y_Vj is 2-dim, and for j != k with
V_j cap V_k = the point p_jk = l_j cap l_k (distance 2), W_j cap W_k = Y_Vj cap Y_Vk =: q_jk (1-dim, canonical).
If V_j cap V_k is a line (distance 1), Y_Vj cap Y_Vk is 2-dim and only constrains W_j cap W_k to lie in it.
Necessary condition used here (exact search over subsets J): the canonical points q_jk (pairs at distance 2 with
l_j != l_k) together span at most 3 dims, and for each j the q_jk (k in J) span at most 2 dims; for l_j = l_k the
2-dim Y_Vj cap Y_Vk must lie in Y_U too.  Reports the largest J passing, per sampled other frame.
"""
import json, sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fermat import P, rank_mod, dec

def nullspace_mod(M):
    """basis of the right null space of M (list of rows) mod P"""
    M = [list(r) for r in M]; R = len(M); C = len(M[0]); piv = []; r = 0
    for c in range(C):
        pv = next((i for i in range(r, R) if M[i][c] % P), None)
        if pv is None: continue
        M[r], M[pv] = M[pv], M[r]; inv = pow(M[r][c], P - 2, P); M[r] = [x * inv % P for x in M[r]]
        for i in range(R):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % P for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
    free = [c for c in range(C) if c not in piv]
    basis = []
    for fc in free:
        v = [0] * C; v[fc] = 1
        for i, pc in enumerate(piv): v[pc] = (-M[i][fc]) % P
        basis.append(v)
    return basis

def intersect(A, B):
    """A, B: lists of row vectors (bases); returns a basis of span(A) cap span(B)"""
    # x A = y B  <=>  [A; -B]^T (x, y) = 0
    rows = A + [[(-x) % P for x in b] for b in B]
    MT = [list(col) for col in zip(*rows)]            # 6 x (a + b)
    ns = nullspace_mod(MT)
    out = []
    for v in ns:
        x = v[:len(A)]
        w = [sum(x[i] * A[i][c] for i in range(len(A))) % P for c in range(6)]
        out.append(w)
    # reduce to independent
    res = []
    for w in out:
        if rank_mod(res + [w]) > len(res): res.append(w)
    return res

if __name__ == "__main__":
    data = json.load(open("frames891.json")); frames = [set(S) for S in data["frames_D"]]
    Dm = np.load("dist891.npy")
    fe = json.load(open("fermat405.json")); lift = {int(k): v for k, v in fe["lift_modp"].items()}
    isF = set(lift)
    others = [i for i in range(891) if i not in isF]
    sample = others[:: max(1, len(others) // int(sys.argv[1]))] if len(sys.argv) > 1 else others[:12]
    best_overall = 0
    for U in sample:
        nbr = [V for V in range(891) if Dm[U, V] == 1 and V in isF]
        lines = {V: frozenset(frames[U] & frames[V]) for V in nbr}
        # pair data
        q = {}; ok_pair = {}
        for a, b in itertools.combinations(nbr, 2):
            I = intersect(lift[a], lift[b])
            q[(a, b)] = I
        best = 0; bestJ = None
        n = len(nbr)
        for mask in range(1, 1 << n):
            J = [nbr[i] for i in range(n) if mask >> i & 1]
            if len(J) <= best: continue
            pts = []; good = True; per = {j: [] for j in J}; extra = []
            for a, b in itertools.combinations(J, 2):
                key = (a, b) if (a, b) in q else (b, a)
                I = q[key]
                if lines[a] == lines[b]:
                    if len(I) != 2: good = False; break
                    extra += I; per[a] += I; per[b] += I
                elif len(I) == 1:
                    pts += I; per[a] += I; per[b] += I
                # distance-1 pairs with different lines: no canonical constraint here (only necessary conditions used)
            if not good: continue
            if pts + extra and rank_mod(pts + extra) > 3: continue
            if any(per[j] and rank_mod(per[j]) > 2 for j in J): continue
            best = len(J); bestJ = J
        dhist = {}
        for a, b in itertools.combinations(nbr, 2):
            d = int(Dm[a, b]); dhist[d] = dhist.get(d, 0) + 1
        same_line = sum(1 for a, b in itertools.combinations(nbr, 2) if lines[a] == lines[b])
        print(f"other frame {U}: 15 Fermat nbrs; pairs by distance {dhist}; pairs on the same line of U {same_line};"
              f" largest J passing the necessary test: {best}", flush=True)
        best_overall = max(best_overall, best)
    print("max over sample:", best_overall)
