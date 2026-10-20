"""ts.py -- w5-twoshared: common code for the shape with TWO shared moves at layer 2 and ONE shared layer-3 move
(W in 3W+1 for W = 2: s_0, s_1, u_1, u_2, v, t_0, t_1).  Register form (GAME.md s5); engine letters (Clifford
monomials) by default.  Notation: s_j = D b_j x_j; u_k = D(a_k0 s_0 + a_k1 s_1 + iota_k); v = D(g_1 u_1 + g_2 u_2 + ...);
t_e = D(mu_e v + psi_e).  p_kj = D a_kj D b_j (layer 2), q_kj = D g_k p_kj (layer 3), y_e = D mu_e."""
import sys, os, random, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
import numpy as np
from asym import G, Prog, DTest, full2, rand_engine, plain
from fam import prog_letter
from famM import midpaths

n = 16
inv = lambda M: np.conj(M).T
def lam(M):
    return int(round(np.log2((np.abs(M) > 1e-9).sum(axis=1).max())))
def is_kernel(M):
    """M in B'K (widened frame of K)"""
    return G.frame_key(M) == G.frm(G.KID)
def same_frame(M1, M2):
    return G.frame_key(M1) == G.frame_key(M2)
def is_free(M):
    return G.is_monomial(M)
def prop(M1, M2):
    """M1 = zeta M2 for a scalar zeta?"""
    i = np.unravel_index(np.argmax(np.abs(M2)), M2.shape)
    z = M1[i] / M2[i]
    return abs(abs(z) - 1) < 1e-9 and np.allclose(M1, z * M2, atol=1e-9)

def route(rng, tries=50):
    """b, a (progress letters), p = D a D b at layer 2, and the geodesic 2-paths (g, mu) from p to K"""
    for _ in range(tries):
        b = prog_letter(rng, [G.ID]); h = G.intern(G.mul(G.D1, b))
        a = prog_letter(rng, [h]); p = G.mul(G.D1, a, G.MATS[h])
        if lam(p) != 2: continue
        paths = midpaths(rng, G.intern(p), nsamp=300)
        if paths: return b, a, p, paths
    return None
