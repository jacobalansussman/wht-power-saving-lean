"""ctl_annih.py -- ms-xw3.  Positive control for the EXCESS accounting of annih2/3/4: append two synthetic columns
N + Z and Z' - N (N a genuine 3-layer column, Z, Z' in L2) to the plain-chain column set; their sum lies in L2 while
neither does, so the test must report EXCESS 1 (and 0 without them)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from xw3lib import Game, menu as MENU
from annih2 import build


def excess(A3, A2):
    U2, s2, _ = np.linalg.svd(A2, full_matrices=False); Q2 = U2[:, s2 > s2[0] * 1e-9]
    A3p = A3 - Q2 @ (Q2.conj().T @ A3)
    sA = np.linalg.svd(A3, compute_uv=False); rA = int((sA > sA[0] * 1e-9).sum())
    sP = np.linalg.svd(A3p, compute_uv=False); rP = int((sP > sP[0] * 1e-9).sum())
    cn = np.linalg.norm(A3p, axis=0) / np.maximum(np.linalg.norm(A3, axis=0), 1e-300)
    triv = np.flatnonzero(cn < 1e-8); rT = 0
    if triv.size:
        sT = np.linalg.svd(A3[:, triv], compute_uv=False); rT = int((sT > sT[0] * 1e-9).sum())
    return rA - rP, rT, rA - rP - rT


G = Game([3, 4], W=1, S=1, seed=5)
M = MENU("xq")
A3, A2 = build(G, M, M, M, "1;1;2", 4, "1")
print("plain:", excess(A3, A2))
rng = np.random.default_rng(1)
N = rng.standard_normal(A3.shape[0]) + 1j * rng.standard_normal(A3.shape[0])   # outside span A3 and L2
Z = A2 @ rng.standard_normal(A2.shape[1]); Z2 = A2 @ rng.standard_normal(A2.shape[1])
print("with planted pair (N+Z, Z2-N), N generic:", excess(np.concatenate([A3, (N + Z)[:, None], (Z2 - N)[:, None]], 1), A2))
