"""v-ms-fano: two counterexamples to statements in ms-fano/PROOF.md Lemma P, checked with explicit complex matrices
at blocks m = 3 and 4 (array of 2^m entries; the identities are group-algebra identities, so they hold at every f).
(A) frustration with TWO pieces on two different frames, via an order-3 address relabelling sigma (a free letter):
    G1 = x1 + x2, G2 = x1 - sigma x2.  Z2 = K (1 - sigma + sigma^2)/2, Z1 = Z2 sigma give Z1 G1 + Z2 G2 = K x1.
    (PROOF.md (d) and the claim: 'the smallest instance spread over frames is exactly ... a base plus the three points
    of a line' -- true only when the ratio group is an elementary abelian 2-group.)
(B) Lemma P(b)'s 'if' direction fails when X_j is not invertible: G1 = D(1+Z)(x1+x2), G2 = D(1+Z)(x1-x2) are
    'frustrated' (same h, opposite signs) but no P1, P2 give P1 G1 + P2 G2 = K x1: the x1-coefficient (P1+P2) D (1+Z)
    is a singular matrix while K is invertible."""
import numpy as np, itertools
def run(m):
    n = 2 ** m
    def X(z):
        A = np.zeros((n, n), complex)
        for x in range(n): A[x ^ z, x] = 1
        return A
    def Zp(c): return np.diag([(-1) ** bin(c & x).count('1') for x in range(n)]).astype(complex)
    def D(z): return (1 + 1j) / 2 * np.eye(n) + (1 - 1j) / 2 * X(z)
    K = np.eye(n, dtype=complex)
    for t in range(m): K = K @ D(1 << t)
    # sigma: cyclic shift of the address bits 0 -> 1 -> 2 -> 0 (bits >= 3 fixed), a permutation matrix of order 3
    def perm(x): b = [(x >> t) & 1 for t in range(m)]; b[0], b[1], b[2] = b[2], b[0], b[1]; return sum(v << t for t, v in enumerate(b))
    S = np.zeros((n, n), complex)
    for x in range(n): S[perm(x), x] = 1
    I = np.eye(n)
    assert np.allclose(np.linalg.matrix_power(S, 3), I) and not np.allclose(S, I)
    Z2 = K @ (I - S + S @ S) / 2; Z1 = Z2 @ S
    # x1-coefficient and x2-coefficient of Z1 (x1 + x2) + Z2 (x1 - S x2)
    c1 = Z1 + Z2; c2 = Z1 - Z2 @ S
    print(f"m={m} (A) x1-coeff == K: {np.allclose(c1, K)}; x2-coeff == 0: {np.allclose(c2, 0)}; "
          f"max denominator in (1 - s + s^2)/2: 2")
    # frames: sigma conjugates D(e1) to D(e2): the x2 term of D(e1)(x1 - sigma x2), namely D(e1) sigma = sigma D(e?) ...
    for t in range(m):
        Dt = D(1 << t); conj = np.linalg.inv(S) @ Dt @ S
        tgt = [u for u in range(m) if np.allclose(conj, D(1 << u))]
        print(f"   sigma^-1 D(e{t}) sigma = D(e{tgt[0] if tgt else '?'})")
    # (B)
    Zc = Zp(1); Dm = D(1)
    Xb = Dm @ (I + Zc)
    print(f"m={m} (B) rank of X = D(1+Z_e1): {np.linalg.matrix_rank(Xb)} of {n}  -> (P1+P2) X = K impossible")
for m in (3, 4): run(m)
