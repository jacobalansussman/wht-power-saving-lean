# Evaluate a commutative-slice transfer triple (A~, B~, Gamma) at EVERY character of the slice group A (not only the
# order-2 ones): chi(t_z) = i^psi(z), psi in the Z_4-span of z -> [z.y].  Program output = I + A C (I - G C)^-1 B,
# C = diag(chi(t_{z_i}) - 1), G = Gamma/2 (Lemma 1).  Target chi(K) I, K = t_e1 ... t_em.
import sys, itertools
import numpy as np
d = np.load(sys.argv[1])
A, B, Gam = d["A"], d["B"], d["Gam"]; dirs = [int(z) for z in d["dirs"]]; W = int(d["W"]); m = int(d["m"])
M = len(dirs)
dot = lambda z, y: bin(z & y).count("1") & 1
ys = list(range(1 << m))
# all distinct psi: Z_4 combinations of l_y
psis = set()
for ns in itertools.product(range(4), repeat=len(ys)):
    psi = tuple(sum(n * dot(z, y) for n, y in zip(ns, ys)) % 4 for z in range(1 << m))
    psis.add(psi)
psis = sorted(psis)
print("characters:", len(psis))
G = Gam  # in the solver's convention X = A_D (I + G_D)^-1 B_D at order 2, i.e. c = -2 and Gamma = G
worst2 = worst4 = 0.0
for psi in psis:
    w = np.array([1j ** psi[z] for z in dirs])
    c = w - 1
    C = np.diag(c)
    # order-2 normalisation: with c = -2 on D, A C (I - (G/2)... ) ; use the general form with Gamma = G:
    out = np.eye(W) + A @ C @ np.linalg.solve(np.eye(M) - (G / 2) @ C, B)
    target = (1j ** (sum(psi[1 << j] for j in range(m)) % 4)) * np.eye(W)
    err = np.abs(out - target).max()
    if all(p % 2 == 0 for p in psi): worst2 = max(worst2, err)
    else: worst4 = max(worst4, err)
print("max error on order-2 characters %.3e ; on the other characters %.3e" % (worst2, worst4))
