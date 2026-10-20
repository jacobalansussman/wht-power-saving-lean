# c11_staircase.py -- an explicit 5-in-14 scheme in the singular subgame's V_56 relaxation (the D-game), built by the
# "staircase": the free parts alpha_t tile C^140; the 7-parts S_t of moves 1..8 lie in the 7-type target T7, move 9
# finishes T7, moves 10..14 complete the tiling. Then the game is simulated with its own update rules
#   N <- N + Alg (E_h - I) Y,  Y = D + P,  D in N, P in Alg^W   (rows: R7 <- R7 + phi(V7), R21 <- R21 + phi(V))
# and the target (E_K - I) e_l in N is checked for every l. Controls: the same simulation must FAIL with 13 moves.
import numpy as np
dat = np.load('dgame.npz'); G = dat['G']; F = dat['F']; Eh = dat['Eh']
rng = np.random.default_rng(2026)
W = 5; n = 28 * W
def rk(M, tol=1e-8):
    if M.size == 0: return 0
    s = np.linalg.svd(M, compute_uv=False); return int((s > tol * max(1, s[0])).sum())
def rowbasis(M, tol=1e-9):
    u, s, vt = np.linalg.svd(M, full_matrices=False); return vt[s > tol * max(1, s[0])]
P = -G / 2                                   # projector onto the (-1)-eigenspace of E_h (G = E_h - I = -2P)
assert np.allclose(P @ P, P, atol=1e-10)
VP = rowbasis(P); VP7 = rowbasis(P[:7]); VF7 = rowbasis(F[:7]); VF21 = rowbasis(F[7:])
print('dims: V^P', VP.shape[0], ' V^P_7', VP7.shape[0], ' V^F_7', VF7.shape[0], ' V^F_21', VF21.shape[0])
def emb(v, l):
    out = np.zeros(n); out[28 * l:28 * l + 28] = v; return out
T7 = np.array([emb(v, l) for l in range(W) for v in VF7]); T21 = np.array([emb(v, l) for l in range(W) for v in VF21])
# realise a block-diagonal P = diag(p7, p21) with v P = mu(v) on V^P_7 (mu given on the basis rows), random elsewhere
def realise(mu_rows):
    """mu_rows: 4 x 28 targets for the 4 basis rows of VP7. Returns 28 x 28 block-diagonal matrix."""
    V7, V21 = VP7[:, :7], VP7[:, 7:]                 # projections (each rank 4)
    p7 = np.linalg.lstsq(V7, mu_rows[:, :7], rcond=None)[0]
    p21 = np.linalg.lstsq(V21, mu_rows[:, 7:], rcond=None)[0]
    # add random parts acting on complements (do not change values on V^P_7)
    N7 = np.eye(7) - np.linalg.pinv(V7) @ V7; N21 = np.eye(21) - np.linalg.pinv(V21) @ V21
    p7 = p7 + N7 @ rng.standard_normal((7, 7)); p21 = p21 + N21 @ rng.standard_normal((21, 21))
    M = np.zeros((28, 28)); M[:7, :7] = p7; M[7:, 7:] = p21
    assert np.allclose(VP7 @ M, mu_rows, atol=1e-8)
    return M
def build(Mmoves, staircase=True):
    Pt = []
    # target coordinates: T7 basis rows; S_t for t<=8: 4 independent combinations of T7 rows
    Tco = rng.standard_normal((35, 35))          # random basis change of T7
    Tb = Tco @ T7
    for t in range(Mmoves):
        if staircase and t < 8:
            S = Tb[4 * t:4 * t + 4]                   # 4 x 140 rows in T7
        elif staircase and t == 8:
            S = np.vstack([Tb[32:35], rng.standard_normal((1, n))])
        else:
            S = rng.standard_normal((4, n))
        Pt.append([realise(S[:, 28 * l:28 * l + 28]) for l in range(W)])
    return Pt
def alpha(Pl, V):          # rows v -> (v P_1, ..., v P_W)
    return np.hstack([V @ Pl[l] for l in range(W)])
def intersect(Arows, Brows, tol=1e-8):
    """basis rows of span(Arows) cap span(Brows)"""
    if Arows.shape[0] == 0 or Brows.shape[0] == 0: return np.zeros((0, Arows.shape[1]))
    M = np.vstack([Arows, -Brows]).T
    u, s, vt = np.linalg.svd(M)
    r = int((s > tol * max(1, s[0])).sum())
    Nul = vt[r:].T                                   # (a; b) with a A - b B = 0
    X = Nul[:Arows.shape[0]].T @ Arows
    return rowbasis(X) if X.shape[0] else X
def simulate(Pt, use_D=True):
    R7 = np.zeros((0, n)); R21 = np.zeros((0, n)); Fsum = np.zeros((0, n)); Tprev = np.zeros((0, n))
    gains = []
    for t in range(len(Pt)):
        A = alpha(Pt[t], VP)
        before = rk(np.vstack([R7, R21])) if R7.shape[0] + R21.shape[0] else 0
        Fall = np.vstack([Fsum, A])
        Tcur = intersect(T7, Fall)
        # new target directions: complement of Tprev inside Tcur
        if Tprev.shape[0]:
            Q = Tcur - (Tcur @ np.linalg.pinv(Tprev)) @ Tprev       # remove Tprev components
            new = rowbasis(Q) if Q.shape[0] else Q
        else:
            new = Tcur
        D = np.zeros((28, n))
        if use_D and new.shape[0] and Fsum.shape[0]:
            X = np.linalg.lstsq(Fall.T, new.T, rcond=None)[0]
            lead = X[Fsum.shape[0]:]                               # 10 x k
            lower = Fsum.T @ X[:Fsum.shape[0]]                     # 140 x k
            vco = lead.T @ VP                                      # k x 28 (in V^P)
            c7 = np.linalg.lstsq(VP7.T, vco.T, rcond=None)[0]      # 4 x k
            badlead = np.linalg.norm(VP7.T @ c7 - vco.T)
            gam = np.linalg.lstsq(c7.T, lower.T, rcond=None)[0].T  # 140 x 4
            D[7:] = np.linalg.lstsq(VP7[:, 7:], gam.T, rcond=None)[0]
        Ymat = D + np.hstack(Pt[t])
        R7 = np.vstack([R7, VP7 @ Ymat]); R21 = np.vstack([R21, VP @ Ymat]); Fsum = Fall; Tprev = Tcur
        gains.append(rk(np.vstack([R7, R21])) - before)
    ok7 = rk(np.vstack([R7, T7])) == rk(R7); ok21 = rk(np.vstack([R21, T21])) == rk(R21)
    return gains, rk(R7), rk(R21), ok7, ok21
for trial in range(3):
    Pt = build(14)
    gains, d7, d21, ok7, ok21 = simulate(Pt)
    print(f'14 moves (staircase), trial {trial}: hook gains per move {gains}; dim R7 {d7}, dim R21 {d21}; T7 in R7: {ok7}; T21 in R21: {ok21}')
# controls: 13 moves cannot reach (hook); generic (non-staircase) 14 moves should fail T7
Pt = build(13); gains, d7, d21, ok7, ok21 = simulate(Pt)
print(f'control 13 moves: gains {gains}, T7 in R7 {ok7}, T21 in R21 {ok21}')
Pt = build(14, staircase=False); gains, d7, d21, ok7, ok21 = simulate(Pt)
print(f'control 14 generic moves: gains {gains}, T7 in R7 {ok7}, T21 in R21 {ok21}')
Pt = build(14); gains, d7, d21, ok7, ok21 = simulate(Pt, use_D=False)
print(f'control staircase without corrections D: T7 in R7 {ok7}')
