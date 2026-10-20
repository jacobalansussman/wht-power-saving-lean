"""geo9.py (v-g-9hyp-model) -- independent re-implementation: E8 lines (hub Hamming numbering), Gram, B0 column space
over F_2 (no orthonormal-label search), puncture sets.  Written from the hub e8units docstring, not from e8geo.py."""
import itertools
import numpy as np

def lines():
    gen = [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[0,1,0,1,0,1,0,1]]
    code = set()
    for c in itertools.product([0,1], repeat=4):
        w = tuple(sum(c[i]*gen[i][k] for i in range(4)) % 2 for k in range(8))
        code.add(w)
    words = sorted(w for w in code if sum(w) == 4)
    assert len(words) == 14
    out = [tuple(2 if k == a else 0 for k in range(8)) for a in range(8)]
    for w in words:
        sup = [k for k in range(8) if w[k]]
        for sg in itertools.product([1,-1], repeat=3):
            x = [0]*8; x[sup[0]] = 1
            for i in range(3): x[sup[i+1]] = sg[i]
            out.append(tuple(x))
    return out

LINES = lines()
LA = np.array(LINES, dtype=np.int64)
GR = (LA @ LA.T) // 2
assert (np.diag(GR) == 2).all() and len(LINES) == 120
# E8 check: all pairwise |G| <= 1 off the diagonal, 63 orthogonal per line
assert all(sorted(set(abs(GR[i, j]) for j in range(120) if j != i)) == [0, 1] for i in range(120))
assert all(sum(1 for j in range(120) if GR[i, j] == 0) == 63 for i in range(120))
P_S2 = [39, 52, 57, 71, 76, 78, 103, 109, 114]
P_WORD = [8, 17, 39, 51, 53, 58, 63, 72, 114]

def f2_colspace_basis(M):
    """basis (list of python-int bitmasks over rows) of the F_2 column space of 0/1 matrix M (rows x cols)."""
    n, m = M.shape
    cols = [sum(int(M[i, j]) << i for i in range(n)) for j in range(m)]
    basis = {}
    for x in cols:
        y = x
        while y:
            h = y.bit_length() - 1
            if h in basis: y ^= basis[h]
            else: basis[h] = y; break
    return list(basis.values())

def unit(name):
    P = {'s2': P_S2, 'word': P_WORD}[name]
    live = [i for i in range(120) if i not in P]
    B0 = ((1 + GR) % 2).astype(np.int64)
    bas_full = f2_colspace_basis(B0)                      # col space of B0 on all 120 lines
    assert len(bas_full) == 9
    # restrict to live rows, as bitmasks over live positions 0..110
    def restrict(x): return sum(((x >> T) & 1) << t for t, T in enumerate(live))
    W0b = [restrict(x) for x in bas_full]
    return P, live, W0b

def f2_span_all(basis):
    out = [0]
    for b in basis: out = out + [x ^ b for x in out]
    return out
