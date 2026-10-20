"""c4: DFT and polynomial multiplication have the same log-saving exponent (iso-lower, 14 Oct 2026).

Both directions, with O(n) overhead and a constant-factor change of length:
 (=>) Bluestein: with w = exp(-2 pi i/n), jk = (j^2 + k^2 - (k-j)^2)/2, so
        (F_n x)_k = w^{k^2/2} * sum_j (x_j w^{j^2/2}) * w^{-(k-j)^2/2},
      i.e. 2n scalings + ONE product of a variable polynomial (degree < n) with a FIXED chirp polynomial
      (degree < 2n), reading n middle coefficients. With the chirp fixed the product is a linear map; an
      arithmetic circuit for general polynomial products linearises at <= 3x size (check c1 (b)).
 (<=) Convolution theorem: product of two polynomials of degree < n = 3 DFTs of length >= 2n - 1
      + 2n - 1 pointwise products.
So L(DFT, length ~n) <= 3 PM(~2n) + O(n) and PM(n) <= 3 L(DFT, length ~2n) + O(n): any bound
N (log N)^(1-s) transfers both ways with the same s (lengths change by a constant factor, which does not
move s). Checked numerically for several lengths, including primes and the length 30 = 2*3*5.
"""
import sys
import numpy as np

rng = np.random.default_rng(7)

def dft(x):
    n = len(x)
    k = np.arange(n)
    return np.exp(-2j * np.pi * np.outer(k, k) / n) @ x

def polymul_naive(a, b):
    return np.convolve(a, b)

def bluestein(x):
    n = len(x)
    j = np.arange(n)
    w_half = np.exp(-1j * np.pi * j ** 2 / n)                 # w^{j^2/2}
    a = x * w_half                                            # variable polynomial, degree < n
    t = np.arange(-(n - 1), n)                                # chirp indices -(n-1)..(n-1)
    chirp = np.exp(1j * np.pi * t ** 2 / n)                   # w^{-t^2/2}: FIXED polynomial, degree < 2n-1
    c = polymul_naive(a, chirp)                               # one polynomial product
    return w_half * c[n - 1:2 * n - 1]                        # middle coefficients, then scale

def polymul_via_dft(a, b):
    L = len(a) + len(b) - 1
    A = dft(np.concatenate([a, np.zeros(L - len(a))]))
    B = dft(np.concatenate([b, np.zeros(L - len(b))]))
    k = np.arange(L)
    inv = np.exp(2j * np.pi * np.outer(k, k) / L) / L
    return inv @ (A * B)

for n in (2, 3, 5, 7, 8, 13, 16, 30, 31):
    x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    e1 = np.abs(bluestein(x) - dft(x)).max()
    a = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    e2 = np.abs(polymul_via_dft(a, b) - polymul_naive(a, b)).max()
    print(f"n={n:3d}: |Bluestein - DFT| = {e1:.1e}   |DFT-based product - naive product| = {e2:.1e}")
