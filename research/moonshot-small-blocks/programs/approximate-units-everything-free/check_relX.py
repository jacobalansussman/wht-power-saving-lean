"""v-w3-border: does EXACTIFY's Theorem X (relative version) keep the chain-map commutation a1 D1 = D1 a1?
Toy border family (one paid symbol X, target D1 c):
    Phi(eps; X) = X c + eps^-1 (a1 X c' - X a1 c')          (3 moves of kind 1: X c, X c', X a1 c')
* evaluated at X = D1: the eps^-1 part is 0 (a1 commutes with D1)  -> a border (here even exact-at-every-eps) family;
* STRICTLY (X a free symbol, w3-border's A = R * C<X>) the eps^-1 part is NOT zero -> by w3-border's section 1 it is
  'relation-using' (F2 case (b): the formal output explodes along the relation a1 X - X a1);
* Theorem X's recursion: the sub-call for D1 returns Y(eps) = D1 + eps Y1 + eps^2 Y2 + ...; it is computed on the f
  bit-0 coordinates and applied to every bit-1 fibre, so Y_k = E'_k (x) I.  Then the order-0 coefficient of
  Phi(eps; Y(eps)) is D1 c exactly (relation kept);  with a generic (non-fibrewise) error it is wrong.
Address: column j bits (b0,b1) at positions 2j, 2j+1.  Dense float64 matrices, f = 2..5 (N <= 1024)."""
import numpy as np
rng = np.random.default_rng(11)
a, b = (1 + 1j) / 2, (1 - 1j) / 2

def perm_a1(f):
    x = np.arange(4 ** f); y = x.copy()
    for j in range(1, f):
        y ^= ((x >> (2 * (j - 1) + 1)) & 1) << (2 * j + 1)
    P = np.zeros((4 ** f, 4 ** f)); P[y, x] = 1; return P

def bit0_embed(Ep, f):
    """E' acting on the f bit-0 coordinates (index sum_j b0_j 2^j), identity on bit-1 coordinates"""
    N = 4 ** f; x = np.arange(N)
    s0 = sum(((x >> (2 * j)) & 1) << j for j in range(f))
    s1 = sum(((x >> (2 * j + 1)) & 1) << j for j in range(f))
    M = Ep[np.ix_(s0, s0)] * (s1[:, None] == s1[None, :])
    return M

def D1(f):
    d = np.array([[a, b], [b, a]])
    T = np.array([[1.0 + 0j]])
    for _ in range(f): T = np.kron(d, T)
    return bit0_embed(T, f)

for f in range(2, 6):
    N = 4 ** f; P = perm_a1(f); D = D1(f)
    c = rng.standard_normal(N) + 1j * rng.standard_normal(N); cp = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    # strict formal check: replace X by a generic matrix (a representation of the free product R * C<X>)
    Xg = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    strict = np.linalg.norm(P @ Xg @ cp - Xg @ P @ cp) / np.linalg.norm(cp)
    ev = np.linalg.norm(P @ D @ cp - D @ P @ cp)
    # Theorem X: Y(eps) = D + eps Y1 + eps^2 Y2; order-0 coefficient of Phi(eps;Y) = D c + (a1 Y1 - Y1 a1) c'
    E1f = bit0_embed(rng.standard_normal((2 ** f, 2 ** f)) + 1j * rng.standard_normal((2 ** f, 2 ** f)), f)
    E1g = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    o0_fib = D @ c + (P @ E1f - E1f @ P) @ cp
    o0_gen = D @ c + (P @ E1g - E1g @ P) @ cp
    # also a depth-2 style check: sub-call itself a Laurent series with negative-free window, fibrewise
    print(f"f={f}: |a1 D1 - D1 a1| c' = {ev:.1e};  strict eps^-1 part (generic X) = {strict:.2f} (nonzero: relation-using);"
          f"  order-0 error, fibrewise sub-call error = {np.linalg.norm(o0_fib - D @ c):.1e};"
          f"  generic error = {np.linalg.norm(o0_gen - D @ c):.1f};  [E1f, a1] = {np.linalg.norm(P @ E1f - E1f @ P):.1e}")
