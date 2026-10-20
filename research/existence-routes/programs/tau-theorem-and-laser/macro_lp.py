# macro_lp.py (meta/ex-laser, 14 Oct 2026).  Route (a) of PROOF.md section 3: the engine's EXISTING pricing.
# The published engine prices a pure-target word at block B with W arrays and grouped paid blocks of ranks r < B by
# the moment root: the exponent tau with  sum_r n_r r^tau = W B^tau  (stages.py root(); 1 - tau = saving).
# A mixed unit (W3 kernels + W2 rank-2 targets in M moves at block 3) can be packaged inside a "macro-word" at block B:
# one B-bit kernel demand per array, split freely, fed to the mixed unit at any column c (3c, 2c -> M demands of c),
# and every demand of size r < B left over becomes a paid block of rank r.  For fixed tau the LP
#     min  sum_r (mass of leftover rank-r blocks) * r^tau     from one B-demand
# is the best macro-word at block B (fractional arrays; integral words approximate it with many arrays).
# tau_B = min tau with LP(tau) <= B^tau.  Prediction: tau_B > tau_c for every B (the top split loses Phi), and
# tau_B -> tau_c like 1/log B.
import sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
from math import log

def lp(B, W3, W2, M, tau):
    var = []
    for k in range(2, B + 1):
        for a in range(1, k // 2 + 1): var.append(('s', k, a))
    for s in range(1, B // 3 + 1): var.append(('u', s))
    for r in range(1, B): var.append(('p', r))          # leave a size-r demand as a paid block of rank r
    idx = {v: i for i, v in enumerate(var)}
    A = lil_matrix((B + 1, len(var))); b = np.zeros(B + 1); c = np.zeros(len(var))
    for v, i in idx.items():
        if v[0] == 's':
            _, k, a = v
            A[k, i] += 1; A[a, i] -= 1; A[k - a, i] -= 1
        elif v[0] == 'u':
            s = v[1]
            A[3 * s, i] += W3; A[2 * s, i] += W2; A[s, i] -= M
        else:
            r = v[1]; A[r, i] += 1; c[i] = r ** tau
    b[B] = 1
    res = linprog(c, A_eq=A[1:].tocsr(), b_eq=b[1:], bounds=(0, None), method='highs')
    assert res.status == 0, res.message
    return res.fun

def tau_B(B, W3, W2, M):
    lo, hi = 0.5, 1.0
    for _ in range(30):
        mid = (lo + hi) / 2
        if lp(B, W3, W2, M, mid) <= B ** mid: hi = mid
        else: lo = mid
    return hi

def tau_c(W3, W2, M):
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if W3 * 3 ** mid + W2 * 2 ** mid > M: hi = mid
        else: lo = mid
    return lo

if __name__ == "__main__":
    Bs = [int(x) for x in sys.argv[1:]] or [3, 6, 12, 24, 48, 96, 192]
    for (W3, W2, M) in ((1, 4, 10), (7, 38, 88), (5, 0, 14)):
        tc = tau_c(W3, W2, M)
        print(f"{W3}+{W2} in {M}: tau_c {tc:.5f} (saving {1-tc:.5f}); conversion saving "
              f"{1 - log((M + W2) / (W3 + W2)) / log(3):.5f}", flush=True)
        for B in Bs:
            t = tau_B(B, W3, W2, M)
            print(f"   macro-word at block B={B}: tau_B {t:.5f}, saving {1-t:.5f}, "
                  f"(tau_B - tau_c)*ln B = {(t - tc) * log(B):.4f}", flush=True)
