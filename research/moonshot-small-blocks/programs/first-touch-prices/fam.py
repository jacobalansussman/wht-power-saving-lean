# fam.py -- w7-mixed3: asymmetric layered skeletons for block-3 mixed units, screened by the f = 2 widened test.
# A skeleton = list of moves (zdir, items) with items (coef, letter, src) as in blk3.py; coefficients are the unknowns.
# Screen: minimise the f = 2 residual (variable projection: the readout is solved exactly per block) over ALL
# coefficients by Levenberg-Marquardt from random starts; judge with the full residual on fresh random coefficients'
# optimum (all three blocks) and, for any low value, on the f = 3 diagonal block.
import sys, os, time, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from blk3 import *
from scipy.optimize import least_squares

def kr(M): return np.kron(M, M)
D2M = {z: kr(D[z]) for z in (1, 2, 4)}
TG2 = {k: kr(v) for k, v in TGT.items()}

class Skel:
    def __init__(self, W, targets, moves, blocks=('D', 'S', 'A')):
        self.W, self.targets, self.moves = W, targets, moves
        self.n = ncoefs((W, targets, moves))
        self.items = []
        for m, (z, its) in enumerate(moves):
            for (c, L, src) in its: self.items.append((m, kr(L), src))
        self.blocks = blocks
        self.Tw = []
        for w in range(W):
            for b in blocks:
                Q = QB[b]
                T = np.zeros((Q.shape[0], 64 * W), complex); T[:, 64 * w:64 * (w + 1)] = Q @ TG2[targets[w]]
                self.Tw.append((b, T / np.linalg.norm(T)))
    def gens(self, c):
        W = self.W
        G = [np.zeros((W, 64, 64), complex) for _ in range(W)]
        for w in range(W): G[w][w] = np.eye(64)
        k = 0
        for m, (z, its) in enumerate(self.moves):
            acc = np.zeros((W, 64, 64), complex)
            for (_, L, src) in its:
                acc += c[k] * np.einsum('ij,ljk->lik', self.items[k][1], G[src]) if False else c[k] * (self.items[k][1] @ G[src])
                k += 1
            G.append(D2M[z] @ acc)
        return G
    def resvec(self, x, return_tot=False):
        c = x[:self.n] + 1j * x[self.n:]
        G = self.gens(c)
        out = []; tot = 0.0
        Sb = {}
        for b in self.blocks:
            Q = QB[b]
            S = np.vstack([np.hstack([Q @ g[l] for l in range(self.W)]) for g in G])
            # orthonormal row basis
            U, s, Vh = np.linalg.svd(S, full_matrices=False)
            r = int((s > 1e-9 * s[0]).sum())
            Sb[b] = Vh[:r]
        for (b, T) in self.Tw:
            V = Sb[b]
            R = T - (T @ V.conj().T) @ V
            out.append(R.ravel()); tot += float(np.sum(np.abs(R) ** 2))
        v = np.concatenate(out)
        if return_tot: return tot
        return np.concatenate([v.real, v.imag])
    def total(self, x): return self.resvec(x, return_tot=True)

def solve(sk, rng, iters=60, x0=None):
    x0 = rng.normal(size=2 * sk.n) if x0 is None else x0
    r = least_squares(sk.resvec, x0, method='lm', max_nfev=iters * (2 * sk.n + 1), xtol=1e-12, ftol=1e-14)
    return r.x, sk.total(r.x)

# ---------------------------------------------------------------- skeleton families ---------------------------------
def _pool():
    P = [I8]
    sig = mono([((x << 1) & 6) | (x >> 2) for x in range(8)], [0] * 8)
    P += [sig, sig @ sig]
    for z in (1, 2, 4): P.append(mono([x ^ z for x in range(8)], [0] * 8))                      # shifts
    for w in (1, 2, 4, 3, 5, 6, 7): P.append(mono(list(range(8)), [2 * (bin(w & x).count('1') & 1) for x in range(8)]))  # Z-phases
    for w in (1, 2, 4, 3, 5, 6, 7): P.append(mono(list(range(8)), [(bin(w & x).count('1') & 1) for x in range(8)]))      # S-phases
    for (a, b) in ((0, 1), (1, 2), (0, 2)):
        P.append(mono(list(range(8)), [2 * ((x >> a) & (x >> b) & 1) for x in range(8)]))      # CZ
        P.append(mono([x ^ (((x >> a) & 1 ^ (x >> b) & 1) * ((1 << a) | (1 << b))) for x in range(8)], [0] * 8))  # swap bits
        P.append(mono([x ^ (((x >> a) & 1) << b) for x in range(8)], [0] * 8))                # CNOT a -> b
        P.append(mono([x ^ (((x >> b) & 1) << a) for x in range(8)], [0] * 8))
    return P
POOL = _pool()

def menu(rng, kind, k):
    """k letters: identity first, then random engine / widened / structured-pool letters"""
    out = [I8]
    for _ in range(k - 1):
        if kind == 'engine': out.append(rand_engine(rng))
        elif kind == 'widened': out.append(rand_widened(rng))
        elif kind == 'struct': out.append(POOL[int(rng.integers(1, len(POOL)))])
        elif kind == 'mix':
            u = rng.random()
            out.append(POOL[int(rng.integers(1, len(POOL)))] if u < 0.5 else rand_engine(rng))
    return out

SIG = mono([((x << 1) & 6) | (x >> 2) for x in range(8)], [0] * 8)     # cyclic shift of the 3 address bits (order 3)

def layered(nK, nR, n1, n2, n3, rng, kind='engine', k1=1, k2=2, k3=2, dirs=(1, 2, 4), inj2=0):
    """dense layered skeleton: n1 layer-1 moves reading every input, n2 layer-2 moves reading every layer-1 move
    (and inj2 letters of every input), n3 layer-3 moves reading every layer-2 and layer-1 move."""
    W = nK + nR
    targets = ['K'] * nK + ['R'] * nR
    moves = []
    L1 = list(range(W, W + n1)); L2 = list(range(W + n1, W + n1 + n2))
    for t in range(n1):
        moves.append((dirs[0], [(1.0, L, l) for l in range(W) for L in menu(rng, kind, k1)]))
    for s in range(n2):
        its = [(1.0, L, t) for t in L1 for L in menu(rng, kind, k2)]
        its += [(1.0, L, l) for l in range(W) for L in menu(rng, kind, inj2)] if inj2 else []
        moves.append((dirs[1], its))
    for q in range(n3):
        its = [(1.0, L, s) for s in L2 for L in menu(rng, kind, k3)] + [(1.0, L, t) for t in L1 for L in menu(rng, kind, 1)]
        moves.append((dirs[2], its))
    return Skel(W, targets, moves)

def krelay(nR, nL1, nL2, rng, kind='engine', k=2):
    """K chain k1 = D1(x0 + R inputs), k2 = D2(k1 + L1 rows), k3 = D4(k2 + L2 rows + L1 rows); R: nL1 layer-1 rows,
    nL2 layer-2 rows reading k1 and all L1 rows.  Moves = 3 + nL1 + nL2."""
    W = 1 + nR; targets = ['K'] + ['R'] * nR
    moves = []
    moves.append((1, [(1.0, L, l) for l in range(W) for L in menu(rng, kind, 1 if l else 1)]))          # k1 = W
    L1 = list(range(W + 1, W + 1 + nL1))
    for t in range(nL1): moves.append((1, [(1.0, L, l) for l in range(W) for L in menu(rng, kind, 1)]))
    k1 = W
    L2 = list(range(W + 1 + nL1, W + 1 + nL1 + nL2))
    for s in range(nL2): moves.append((2, [(1.0, L, t) for t in [k1] + L1 for L in menu(rng, kind, k)]))
    k2i = W + 1 + nL1 + nL2
    moves.append((2, [(1.0, L, t) for t in [k1] + L1 for L in menu(rng, kind, k)]))                     # k2
    moves.append((4, [(1.0, L, t) for t in [k2i] + L2 for L in menu(rng, kind, k)] + [(1.0, L, t) for t in [k1] + L1 for L in menu(rng, kind, 1)]))
    return Skel(W, targets, moves)


def krelayN(nK, nR, nL1, nL2, rng, kind='engine', k=2, kin=1):
    """nK kernel chains k1_i = D1(x_Ki + every R input [kin letters]), k2_i = D2(k1's + L1 rows), k3_i = D4(k2's + L2 rows +
    k1's + L1 rows); R arrays: nL1 layer-1 rows (every input), nL2 layer-2 rows reading every k1 and L1 row.
    Moves = 3 nK + nL1 + nL2."""
    W = nK + nR; targets = ['K'] * nK + ['R'] * nR
    moves = []
    for i in range(nK):
        moves.append((1, [(1.0, I8, i)] + [(1.0, L, l) for l in range(nK, W) for L in menu(rng, kind, kin)]))
    K1 = list(range(W, W + nK))
    L1 = list(range(W + nK, W + nK + nL1))
    for t in range(nL1): moves.append((1, [(1.0, L, l) for l in range(W) for L in menu(rng, kind, 1)]))
    L2 = list(range(W + nK + nL1, W + nK + nL1 + nL2))
    for s_ in range(nL2): moves.append((2, [(1.0, L, t) for t in K1 + L1 for L in menu(rng, kind, k)]))
    K2 = list(range(W + nK + nL1 + nL2, W + 2 * nK + nL1 + nL2))
    for i in range(nK): moves.append((2, [(1.0, L, t) for t in K1 + L1 for L in menu(rng, kind, k)]))
    for i in range(nK):
        moves.append((4, [(1.0, L, t) for t in K2 + L2 for L in menu(rng, kind, k)] + [(1.0, L, t) for t in K1 + L1 for L in menu(rng, kind, 1)]))
    return Skel(W, targets, moves)


def dag(nK, nR, M, rng, kind='struct', k=2, dirs='cyc'):
    """fully dense DAG: move m reads every input and every earlier generator with k letters (identity first);
    directions cycle 1, 2, 4 (dirs='cyc') or are random (dirs='rand')."""
    W = nK + nR; targets = ['K'] * nK + ['R'] * nR
    moves = []
    for m in range(M):
        z = (1, 2, 4)[m % 3] if dirs == 'cyc' else (1, 2, 4)[int(rng.integers(3))]
        moves.append((z, [(1.0, L, src) for src in range(W + m) for L in menu(rng, kind, k)]))
    return Skel(W, targets, moves)


def layeredmix(nK, nR, n1, n2, n3, rng, kind='struct', k2=2, inj2=0, d2=(2, 4), d3=(4, 2)):
    """layered skeleton whose middle-layer moves alternate directions d2 and whose top moves alternate d3, so that
    two-row arrays can be finished in different orientations (D1 D2 or D1 D4 up to letters)"""
    W = nK + nR; targets = ['K'] * nK + ['R'] * nR
    moves = []
    L1 = list(range(W, W + n1)); L2 = list(range(W + n1, W + n1 + n2))
    for t in range(n1): moves.append((1, [(1.0, L, l) for l in range(W) for L in menu(rng, kind, 1)]))
    for s_ in range(n2):
        its = [(1.0, L, t) for t in L1 for L in menu(rng, kind, k2)]
        its += [(1.0, L, l) for l in range(W) for L in menu(rng, kind, inj2)] if inj2 else []
        moves.append((d2[s_ % len(d2)], its))
    for q in range(n3):
        its = [(1.0, L, s_) for s_ in L2 for L in menu(rng, kind, k2)] + [(1.0, L, t) for t in L1 for L in menu(rng, kind, 1)]
        moves.append((d3[q % len(d3)], its))
    return Skel(W, targets, moves)

# ---------------------------------------------------------------- analytic gradient (adjoint) ------------------------
def fgrad(sk, x):
    """F = sum_{w,b} ||T_wb (I - P_b)||^2 (normalised targets) and dF/dx, x = (Re c, Im c)"""
    W = sk.W; n = sk.n
    c = x[:n] + 1j * x[n:]
    G = sk.gens(c)
    ng = len(G)
    Gbar = [np.zeros((W, 64, 64), complex) for _ in range(ng)]
    F = 0.0
    for b in sk.blocks:
        Q = QB[b]; d = Q.shape[0]
        S = np.vstack([np.hstack([Q @ g[l] for l in range(W)]) for g in G])      # (ng*d) x 64W
        U, s, Vh = np.linalg.svd(S, full_matrices=False)
        r = int((s > 1e-9 * s[0]).sum())
        Ur, sr, Vr = U[:, :r], s[:r], Vh[:r]                                        # P = Vr^H Vr
        Spinv = Vr.conj().T @ (Ur.conj().T / sr[:, None])                          # 64W x rows
        Msum = np.zeros((64 * W, 64 * W), complex)
        for (bb, T) in sk.Tw:
            if bb != b: continue
            R = T - (T @ Vr.conj().T) @ Vr
            F += float(np.sum(np.abs(R) ** 2))
            Msum += T.conj().T @ R                                                  # T^H T (I - P)
        # dF = -2 Re tr( (I-P) T^H T S^+ dS ) = 2 Re <Sbar, dS>, Sbar = -( (I-P) T^H T S^+ )^H = -(S^+)^H Msum^H ... careful
        Gam = Msum.conj().T @ Spinv            # ((I-P)^H... ) : (I-P) T^H T  = Msum^H since (I-P) Hermitian
        Sbar = -Gam.conj().T                   # rows x 64W
        for j in range(ng):
            blk = Sbar[j * d:(j + 1) * d]
            for l in range(W):
                Gbar[j][l] += Q.conj().T @ blk[:, 64 * l:64 * (l + 1)]
    grad_c = np.zeros(n, complex)    # stores <Abar, L G_src>
    # backward through moves
    k_end = n
    for m in range(len(sk.moves) - 1, -1, -1):
        z, its = sk.moves[m]
        gi = W + m
        Abar = D2M[z].conj().T @ Gbar[gi]
        k0 = k_end - len(its)
        for kk in range(k0, k_end):
            _, L2, src = sk.items[kk]
            LG = L2 @ G[src]
            grad_c[kk] = np.vdot(Abar.ravel(), LG.ravel())
            Gbar[src] += np.conj(c[kk]) * (L2.conj().T @ Abar)
        k_end = k0
    g = np.concatenate([2 * grad_c.real, -2 * grad_c.imag])
    return F, g

def scaled_x0(sk, rng, scale=1.0):
    """random start with each move's coefficients scaled by 1/sqrt(#items) (keeps generator norms O(1))"""
    re = np.concatenate([rng.normal(size=len(its)) * scale / np.sqrt(len(its)) for (z, its) in sk.moves])
    im = np.concatenate([rng.normal(size=len(its)) * scale / np.sqrt(len(its)) for (z, its) in sk.moves])
    return np.concatenate([re, im])

def solve_bfgs(sk, rng, maxiter=400, x0=None, scaled=False):
    from scipy.optimize import minimize
    if x0 is None: x0 = scaled_x0(sk, rng) if scaled else rng.normal(size=2 * sk.n)
    r = minimize(lambda x: fgrad(sk, x), x0, jac=True, method='L-BFGS-B', options={'maxiter': maxiter, 'gtol': 1e-12, 'ftol': 1e-15})
    return r.x, r.fun

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('fam'); ap.add_argument('--seed', type=int, default=1); ap.add_argument('--trials', type=int, default=3)
    ap.add_argument('--kind', default='engine'); ap.add_argument('--iters', type=int, default=400)
    ap.add_argument('--shape', default='6,5,1'); ap.add_argument('--nK', type=int, default=1); ap.add_argument('--nR', type=int, default=5)
    ap.add_argument('--k2', type=int, default=2); ap.add_argument('--blocks', default='DSA'); ap.add_argument('--inj2', type=int, default=0); ap.add_argument('--kin', type=int, default=1); ap.add_argument('--dirs', default='cyc')
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    for tr in range(a.trials):
        t0 = time.time()
        if a.fam == 'layered':
            n1, n2, n3 = map(int, a.shape.split(','))
            sk = layered(a.nK, a.nR, n1, n2, n3, rng, a.kind, k2=a.k2, inj2=a.inj2)
        elif a.fam == 'krelay':
            nL1, nL2 = map(int, a.shape.split(','))
            sk = krelay(a.nR, nL1, nL2, rng, a.kind, k=a.k2)
        elif a.fam == 'layeredmix':
            n1, n2, n3 = map(int, a.shape.split(','))
            sk = layeredmix(a.nK, a.nR, n1, n2, n3, rng, a.kind, k2=a.k2, inj2=a.inj2)
        elif a.fam == 'dag':
            sk = dag(a.nK, a.nR, int(a.shape), rng, a.kind, k=a.k2, dirs=a.dirs)
        elif a.fam == 'krelayN':
            nL1, nL2 = map(int, a.shape.split(','))
            sk = krelayN(a.nK, a.nR, nL1, nL2, rng, a.kind, k=a.k2, kin=a.kin)
        sk.blocks = tuple(a.blocks)
        sk.Tw = [t for t in sk.Tw if t[0] in sk.blocks]
        x, tot = solve_bfgs(sk, rng, a.iters, scaled=(a.fam == 'dag'))
        # judge on all blocks
        full = Skel(sk.W, sk.targets, sk.moves); full_tot = full.total(x)
        print(json.dumps({'fam': a.fam, 'shape': a.shape, 'kind': a.kind, 'seed': a.seed, 'trial': tr, 'moves': len(sk.moves),
                          'params': 2 * sk.n, 'fit_res': round(tot, 6), 'all_blocks_res': round(full_tot, 6),
                          'secs': round(time.time() - t0, 1)}), flush=True)
        np.save(f'x_{a.fam}_{a.shape.replace(",", "")}_{a.kind}_s{a.seed}_t{tr}_k{a.k2}_i{a.inj2}_{a.nK}{a.nR}.npy', x)

