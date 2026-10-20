"""sp_check.py -- validate sp.py against exact matrices (bm.Game): layers and frame equality for random words"""
import sys, os, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import numpy as np
from asym import G, rand_engine
import sp
# Pauli matrices for (x|z): X^x Z^z
def pauli(v):
    M = np.eye(16, dtype=complex)
    x = v & 15; z = v >> 4
    X = np.zeros((16, 16)); 
    for a in range(16): X[a ^ x, a] = 1
    Zd = np.diag([(-1) ** bin(a & z).count('1') for a in range(16)])
    return X @ Zd
def symp_of(g):
    S = []
    for i in range(8):
        P = pauli(1 << i); Q = g @ P @ np.conj(g).T
        found = None
        for v in range(256):
            R = pauli(v); 
            k = np.argmax(np.abs(R).ravel()); ph = Q.ravel()[k] / R.ravel()[k]
            if abs(abs(ph) - 1) < 1e-9 and np.allclose(Q, ph * R, atol=1e-9): found = v; break
        assert found is not None
        S.append(found)
    return tuple(S)
rng = random.Random(5)
SD = symp_of(G.D1); assert SD == sp.D, (SD, sp.D)
SK = symp_of(G.K); assert sp.lag(SK) == sp.LK
bad = 0
for t in range(60):
    g = np.eye(16, dtype=complex); S = sp.I8
    for k in range(rng.randrange(1, 6)):
        a = rand_engine(rng); g = G.mul(G.D1, a, g)
    S = symp_of(g)
    lam_exact = int(round(np.log2((np.abs(g) > 1e-9).sum(axis=1).max())))
    if lam_exact != sp.lam(S): bad += 1
    # frame equality against another word: compare frame keys vs lagrangians
    h = G.mul(rand_engine(rng), g)
    assert sp.lag(symp_of(h)) == sp.lag(S)
    assert G.frame_key(h) == G.frame_key(g)
print("layer mismatches", bad, "of 60; D, K checked")
# random letter generator: symplectic and stabilises L0
for t in range(200):
    a = sp.letter(rng); assert sp.is_symp(a)
print("letters ok")
