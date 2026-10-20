import itertools, numpy as np, sys
from vlib import *
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
L = lagrangians(); say("Lagrangians:", len(L))
E = [(i, j) for i in range(len(L)) for j in range(i + 1, len(L)) if len(L[i] & L[j]) == 4]
say("edges (meet of dim 2):", len(E))
iZ = L.index(LZ)
d1 = [j for j in range(len(L)) if len(L[j] & LZ) == 4]; say("cells at distance 1 from L_Z:", len(d1))
sK = sym(Kmat); LY = apply_sym(sK, LZ); iY = L.index(LY)
say("L_Y = sigma(K) L_Z meets L_Z in", len(LY & LZ), "element(s) (distance 3 iff 1)")
# Levi lifts: P_A * diag phase * Pauli, free, with K^-1 g K free
def dq(Sm):
    d = []
    for x in range(8):
        v = [(x >> j) & 1 for j in range(3)]
        q = sum(Sm[j][j] * v[j] for j in range(3)) + 2 * sum(Sm[j][k] * v[j] * v[k] for j in range(3) for k in range(j + 1, 3))
        d.append(1j ** (q % 4))
    return np.diag(d)
SYMS = [[[b[0], b[3], b[4]], [b[3], b[1], b[5]], [b[4], b[5], b[2]]] for b in itertools.product([0, 1], repeat=6)]
def order_mod_scalar(g, maxo=50):
    h = g.copy()
    for k in range(1, maxo):
        if np.allclose(h, h[0, 0] * np.eye(8)) and abs(h[0, 0]) > 1e-9: return k, h[0, 0]
        h = h @ g
    return None, None
def levi_lifts(A):
    PA = permm(lambda x: lin(A, x)); res = []
    for Sm in SYMS:
        g = PA @ dq(Sm)
        if is_free(Kinv @ g @ Kmat): res.append(g)
    return res
# Singer: A of order 7
G3 = gl3()
def ordA(A):
    f = lambda x: x; k = 0
    cur = list(range(8))
    for k in range(1, 30):
        cur = [lin(A, y) for y in cur]
        if cur == list(range(8)): return k
o7 = [A for A in G3 if ordA(A) == 7]; o3 = [A for A in G3 if ordA(A) == 3]
say("|GL(3,2)| =", len(G3), " order 7:", len(o7), " order 3:", len(o3))
A7 = o7[0]
lifts = levi_lifts(A7); say("Levi lifts (phase-diag) of one Singer A:", len(lifts))
c = lifts[0]; k, s = order_mod_scalar(c); c = c / s ** (1 / 7)
assert np.allclose(np.linalg.matrix_power(c, 7), np.eye(8)); say("Singer lift c: exact order 7 after scalar; free:", is_free(c), " K^-1 c K free:", is_free(Kinv @ c @ Kmat))
sc = sym(c)
fixc = [i for i in range(len(L)) if apply_sym(sc, L[i]) == L[i]]
fixe = [e for e in E if {L.index(apply_sym(sc, L[e[0]])), L.index(apply_sym(sc, L[e[1]]))} == set(e)]
say("Singer fixes cells:", len(fixc), "(L_Z, L_Y among them:", iZ in fixc, iY in fixc, ") fixes edges:", len(fixe))
# every order-7 element of P0 = Stab(L_Z): linear part on L_Z; invariant 2-subspaces of L_Z
inv2 = 0
for A in o7:
    # A acts on z-part by A^-T for a permutation lift; any order-7 linear map: check invariant hyperplanes
    for H in itertools.combinations(range(1, 8), 2):
        Hs = span(H)
        if len(Hs) == 4 and frozenset(lin(A, x) for x in Hs) == Hs: inv2 += 1
say("invariant 2-dim subspaces summed over the 48 order-7 elements of GL(3,2):", inv2)
# hook calibration
cap = np.linalg.matrix_rank(np.hstack([Nb, rho(dirm(1)) @ Nb]), tol=1e-7) - 280
EN = cap_meet(Nb, rho(Kinv) @ Nb)
say("hook: cap(dir1) =", cap, " dim E_N =", EN.shape[1], " demand =", 280 - EN.shape[1])
for z in range(1, 8):
    assert np.linalg.matrix_rank(np.hstack([Nb, rho(dirm(z)) @ Nb]), tol=1e-7) - 280 == 80
say("cap 80 for every dir z, z = 1..7")
cp = Kinv @ inv(c) @ Kmat; assert is_free(cp)
Rc = [rho(np.linalg.matrix_power(cp, s)) for s in range(7)]
say("c' on N/E_N multiplicities:", mults([quot_trace(R, Nb, EN) for R in Rc], 7))
np.save("v_c.npy", c); np.save("v_EN.npy", EN)
open("out_v1_basic.txt", "w").write("\n".join(out) + "\n")
