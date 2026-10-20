"""f2blind.py -- e3-bpexact: the f = 2 widened test (w3-asym rep2.py) with the move contents left as UNKNOWNS is blind.

Explicit construction [computed]: ONE array, TWO moves D(z1), D(z2), both reading the input with widened coefficients
A_i (block-diagonal on V(x)V = D16 + S120 + A120), readout R_x, R_1, R_2 (block-diagonal).  We solve
      K(x)K = R_x + R_1 (D1 (x) D1) A_1 + R_2 (D2 (x) D2) A_2           on V (x) V  (256 x 256)
exactly by the rank argument (each move's D16<->S120 off-diagonal blocks have rank 8; K(x)K's have rank <= 16), then
check the identity on FRESH random inputs and with rep2.Rep2.residual (row-space test with these contents fixed).
A program with 2 moves for 1 array (block 4 needs 4: plain is optimal for one array) passes, so a pass of the
survivor in this test carries no information, and the test can never fail on a program whose top orbit reads its
inputs with free contents.  (The engine coefficient algebra on S120 and A120 is also full, so restricting the contents
to engine letters at f = 2 does not change this; argued, see NOTE.)"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rep2
import importlib.util as _ilu
_sp = _ilu.spec_from_file_location('vgrp', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'w5lib', 'vlib', 'grp.py'))
_m = _ilu.module_from_spec(_sp); _sp.loader.exec_module(_m); Grp = _m.Grp
G = Grp(4)
Q = rep2.QB; QD, QS, QA = Q['D'], Q['S'], Q['A']
K2 = np.kron(G.K, G.K)
def blocks(M):
    return {(a, b): Q[a] @ M @ Q[b].T for a in 'DSA' for b in 'DSA'}
KB = blocks(K2)
print('off-diagonal blocks of K(x)K: rank DS', np.linalg.matrix_rank(KB['D', 'S']), ' SD', np.linalg.matrix_rank(KB['S', 'D']),
      ' DA, AD, SA, AS norms', [round(float(np.abs(KB[p]).max()), 12) for p in (('D', 'A'), ('A', 'D'), ('S', 'A'), ('A', 'S'))])
zs = [1, 2]
Ds = [np.kron(G.D[z], G.D[z]) for z in zs]
DB = [blocks(M) for M in Ds]
for i, z in enumerate(zs):
    print(f'D({z})(x)D({z}): rank DS', np.linalg.matrix_rank(DB[i]['D', 'S']), ' SD', np.linalg.matrix_rank(DB[i]['S', 'D']))
def split_rank(M, parts):
    U, s, Vh = np.linalg.svd(M); r = int((s > 1e-10 * s[0]).sum())
    out = []; per = int(np.ceil(r / parts))
    for p in range(parts):
        idx = list(range(p * per, min(r, (p + 1) * per)))
        out.append((U[:, idx] * s[idx], Vh[idx]))
    return out
def factor(M):     # M = Lf Rt with full column / row rank
    U, s, Vh = np.linalg.svd(M); r = int((s > 1e-10 * s[0]).sum())
    return U[:, :r] * s[:r], Vh[:r]
RD = []; RS = []; AD = []; AS = []
for i, (Ui, Wi) in enumerate(split_rank(KB['D', 'S'], 2)):          # K_DS = sum_i R_i^D M_DS^(i) A_i^S
    Lf, Rt = factor(DB[i]['D', 'S'])
    k = Ui.shape[1]
    Rd = np.zeros((16, 16), dtype=complex); Rd[:, :] = np.hstack([Ui, np.zeros((16, Lf.shape[1] - k))]) @ np.linalg.pinv(Lf)
    As = np.linalg.pinv(Rt) @ np.vstack([Wi, np.zeros((Rt.shape[0] - k, 120))])
    RD.append(Rd); AS.append(As)
for i, (Ui, Wi) in enumerate(split_rank(KB['S', 'D'], 2)):          # K_SD = sum_i R_i^S M_SD^(i) A_i^D
    Lf, Rt = factor(DB[i]['S', 'D'])
    k = Ui.shape[1]
    Rs = np.hstack([Ui, np.zeros((120, Lf.shape[1] - k))]) @ np.linalg.pinv(Lf)
    Ad = np.linalg.pinv(Rt) @ np.vstack([Wi, np.zeros((Rt.shape[0] - k, 16))])
    RS.append(Rs); AD.append(Ad)
def blockdiag(d, s, a):
    return Q['D'].T @ d @ Q['D'] + Q['S'].T @ s @ Q['S'] + Q['A'].T @ a @ Q['A']
Z = np.zeros((120, 120))
Rm = [blockdiag(RD[i], RS[i], Z) for i in range(2)]
Am = [blockdiag(AD[i], AS[i], Z) for i in range(2)]
moves = [Ds[i] @ Am[i] for i in range(2)]
rest = K2 - sum(Rm[i] @ moves[i] for i in range(2))
RB = blocks(rest)
print('after the two moves, off-diagonal blocks left (must be 0):', max(float(np.abs(RB[p]).max()) for p in RB if p[0] != p[1]))
Rx = blockdiag(RB['D', 'D'], RB['S', 'S'], RB['A', 'A'])
prog = lambda x: Rx @ x + sum(Rm[i] @ (moves[i] @ x) for i in range(2))
rng = np.random.default_rng(12345)
X = rng.standard_normal((256, 5)) + 1j * rng.standard_normal((256, 5))
print('FRESH random inputs: max |program(x) - (K(x)K) x| / |x| = %.2e' % (np.abs(prog(X) - K2 @ X).max() / np.abs(X).max()))
# rep2's own row-space residual with these two contents FIXED (generators: input, move 1, move 2)
T = rep2.Rep2(1, G.K)
def as_terms(M):    # rep2 wants terms (l, M16, c) and applies rho2 itself; feed the 256x256 directly instead
    return M
blocksG = []
for M in [np.eye(256)] + moves:
    blocksG.append({b: Qb @ M for b, Qb in Q.items()})
tot = 0
for b, Qb in Q.items():
    S = np.vstack([bl[b] for bl in blocksG]); Tt = Qb @ K2
    U, s, Vh = np.linalg.svd(S, full_matrices=False); r = int((s > 1e-9 * s[0]).sum()); Vr = Vh[:r]
    R = Tt - (Tt @ Vr.conj().T) @ Vr; res = float(np.sum(np.abs(R) ** 2) / np.sum(np.abs(Tt) ** 2)); tot += res
    print(f'  rep2-style residual, block {b}: {res:.2e} (row rank {r})')
print('TOTAL f = 2 widened residual of the 1-array 2-move program: %.2e  (0 = passes)' % tot)
# control: the same with ONE move only (rank 8 < 16: must fail)
blocks1 = blocksG[:2]; tot1 = 0
for b, Qb in Q.items():
    S = np.vstack([bl[b] for bl in blocks1]); Tt = Qb @ K2
    U, s, Vh = np.linalg.svd(S, full_matrices=False); r = int((s > 1e-9 * s[0]).sum()); Vr = Vh[:r]
    R = Tt - (Tt @ Vr.conj().T) @ Vr; tot1 += float(np.sum(np.abs(R) ** 2) / np.sum(np.abs(Tt) ** 2))
print('control, input + move 1 only: residual %.3f (must be > 0)' % tot1)
