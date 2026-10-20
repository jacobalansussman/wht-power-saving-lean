"""lemD.py -- w7-xwborder.  Numerical check of the two analytic inputs of Theorem 2W (block 2), at f = 2..7:
 (a) a single paid move, in any free frame (spread move C D_z C^-1, C from a menu with chain maps, CS, boundary and
     parity letters), is orthogonal to every free word in the normalised trace: |tr(w^* C D_z C^-1)|/N <= 2^(-f/2);
 (b) the target K = D1 D2 is orthogonal to every one-paid-letter word c D_z c' (Lemma O, already checked by
     v-mw2-xwlayers: out_lemO.txt).
Traces are computed exactly from the monomial structure: tr(w^* X) = sum_i conj(w[s(i), i]) X[s(i), i]."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import xw2lib

MENU = ["1", "C01", "CS", "P", "a0", "a1", "a0i", "a1i", "A", "Ai", "B", "b1", "E01", "Pr0", "Pr1", "lCS", "fX2", "lZ1",
        "a1.fX2.a1i", "A.CS.Ai", "Ai.CS.A", "CS.C01", "S1", "S3"]

def dense_mono(g, n):
    M = np.zeros((n, n), complex)
    M[g.sig, np.arange(n)] = xw2lib.PH[g.ph]
    return M

for f in range(2, 8):
    G = xw2lib.Game([f]); n = 4 ** f
    worst = 0.0; worstK = 0.0
    # paid matrices D_z (dense, small f only up to 6 to keep memory; f=7: 16384^2 too big -> columns sample)
    if f <= 4:
        D = {z: G.apply_D(z, np.eye(n, dtype=complex), f) for z in (1, 2, 3)}
        Kmat = G.K(np.eye(n, dtype=complex), f)
        mon = {w: dense_mono(G.mono(w)[f], n) for w in MENU}
        for z in (1, 2, 3):
            for c in MENU:
                Cm = mon[c]; Dk = Cm @ D[z] @ Cm.conj().T
                for w in MENU:
                    worst = max(worst, abs(np.trace(mon[w].conj().T @ Dk)) / n)
                    worstK = max(worstK, abs(np.trace(Kmat.conj().T @ (mon[w] @ Dk))) / n)
        print("f=%d  max|tau(w* C D C^-1)| = %.4e   max|tau(K* w C D C^-1)| = %.4e   2^(-f/2) = %.4e"
              % (f, worst, worstK, 2 ** (-f / 2)), flush=True)
    else:
        rng = np.random.default_rng(f)
        cols = rng.choice(n, 64, replace=False)
        E = np.zeros((n, 64), complex); E[cols, np.arange(64)] = 1
        for z in (1, 2, 3):
            for c in MENU:
                g = G.mono(c)[f]
                V = G.apply_mono(g.inv(), E); V = G.apply_D(z, V, f); V = G.apply_mono(g, V)   # C D C^-1 e_j
                for w in MENU:
                    h = G.mono(w)[f]
                    # diagonal entries of w^* (C D C^-1) at sampled columns j: conj(w[s(j), j]) * (CDC^-1)[s(j), j]
                    vals = np.conj(xw2lib.PH[h.ph[cols]]) * V[h.sig[cols], np.arange(64)]
                    worst = max(worst, np.max(np.abs(vals)))
        print("f=%d  max sampled |diagonal entry of w* C D C^-1| = %.4e  (bound on |tau|; 2^(-f/2) = %.4e)"
              % (f, worst, 2 ** (-f / 2)), flush=True)
