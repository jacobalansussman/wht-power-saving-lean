"""dbg2.py: (1) does the engine group contain scalars outside Z4 (e.g. e^{i pi/4})?  (2) for b15b s902, does the base
stream's term at A1 give an EXACT z-producer of the leftover L2 (relative element monomial with Z4 entries)?"""
import numpy as np, random, os
os.environ.setdefault('ZTRIES', '300')
import jdesign as J, split as S, endsolve as es, designer as dz
from sp2 import *
E = es.Eng(M, seed=1)
# (1) scalar test: products of D on one qubit with the phase letter S0 = diag(1, i) on bit 0
n = E.n
S0 = np.diag([1j if (x & 1) else 1 for x in range(n)]).astype(complex)
D = E.D
for name, W in (("(S D)^3", np.linalg.matrix_power(S0 @ D, 3)), ("(D S)^3", np.linalg.matrix_power(D @ S0, 3)),
                ("(S D S)^2 ", np.linalg.matrix_power(S0 @ D @ S0, 2))):
    d = np.diag(W); off = np.abs(W - np.diag(d)).max()
    print(name, "offdiag", round(off, 9), "diag values", sorted(set(np.round(d, 6)))[:4])
