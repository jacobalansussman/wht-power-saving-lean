"""rep2.py -- w3-asym: representation tests of the unit identity at f = 2 (V (x) V, 256-dim), WIDENED free group
B' = Z4 wr S16 (monomial matrices with entries in {1,i,-1,-i}).  Under B', V(x)V = D16 + S120 + A120 (diagonal,
symmetric off-diagonal, antisymmetric off-diagonal), each irreducible with multiplicity 1, so the image algebra of C[B'] is
M16 + M120 + M120 and the left C[B']-span of generators u_j = sum_l a_jl x_l, read at f = 2, is: in block i, the row
space of [Q_i rho(a_j0) | ... | Q_i rho(a_j,W-1)] summed over j.  [K]x_l0 is in the span iff for every block i the rows of
[0..Q_i rho(K)..0] lie in that row space.  This is a NECESSARY condition for a unit (engine and widened game alike,
since engine letters are widened letters).  Residual = squared distance of the target rows from the row space,
divided by the target's squared norm, summed over blocks (0 = passes f = 2)."""
import numpy as np

n = 16
def blocks16():
    Q = {}
    D = np.zeros((16, 256)); 
    for x in range(16): D[x, x * 16 + x] = 1
    S = []; A = []
    for x in range(16):
        for y in range(x + 1, 16):
            v = np.zeros(256); v[x * 16 + y] = 1 / np.sqrt(2); v[y * 16 + x] = 1 / np.sqrt(2); S.append(v)
            w = np.zeros(256); w[x * 16 + y] = 1 / np.sqrt(2); w[y * 16 + x] = -1 / np.sqrt(2); A.append(w)
    return {'D': D, 'S': np.array(S), 'A': np.array(A)}
QB = blocks16()

def rho2(M):
    return np.kron(M, M)

class Rep2:
    """gens: list of dicts {(l, M_key): coef} given as list of lists [(l, matrix, coef)]"""
    def __init__(self, W, K):
        self.W = W; self.K2 = rho2(K)
    def gen_blocks(self, terms):
        """terms: list of (l, M, c) -> dict block -> (d_i x 256W) matrix"""
        X = [np.zeros((256, 256), dtype=complex) for _ in range(self.W)]
        for (l, M, c) in terms:
            X[l] += c * rho2(M)
        return {b: np.hstack([Q @ X[l] for l in range(self.W)]) for b, Q in QB.items()}
    def residual(self, gens_terms, l0, tol=1e-9, detail=False):
        out = {}; tot = 0
        blocks = [self.gen_blocks(t) for t in gens_terms]
        for b, Q in QB.items():
            S = np.vstack([bl[b] for bl in blocks])
            T = np.zeros((Q.shape[0], 256 * self.W), dtype=complex); T[:, 256 * l0:256 * (l0 + 1)] = Q @ self.K2
            U, s, Vh = np.linalg.svd(S, full_matrices=False)
            r = int((s > tol * s[0]).sum()) if s.size else 0
            Vr = Vh[:r]
            R = T - (T @ Vr.conj().T) @ Vr
            res = float(np.sum(np.abs(R) ** 2) / np.sum(np.abs(T) ** 2))
            out[b] = (res, r); tot += res
        return (tot, out) if detail else tot
