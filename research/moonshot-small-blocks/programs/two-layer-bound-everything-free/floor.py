"""floor.py -- v-mw2-xwlayers (skeptic).  Independent test of Theorem B (mw2-xwlayers PROOF.md s5) at block 2, W = 2, M = 3.

Theorem B: for EVERY FIXED coefficient vector a, liminf_f E_f(a) >= W - min(t, c), i.e. relative residual
sqrt(E/W) >= sqrt(1/2) = 0.7071 for shapes with one deep move (ST) or one co-deep move with input (SR).
Test: fit a at a small f window (block coordinate descent, exact least squares per block, own code), then evaluate the
SAME a at f = 2..8 on fresh probe vectors.  A fixed-a residual that stays clearly below 0.707 as f grows refutes B.
CHW (all moves read everything; passes the shape conditions) is the comparison: no prediction.

usage: floor.py SHAPE z1,z2,z3 FITFS(e.g. 2,3) sweeps seed [evalmax]
"""
import sys, time
import numpy as np
from mini import Game

MENU = ["1", "CS", "CSi", "A", "P0", "P1", "C01", "S1", "X0", "chi0", "chi1i", "lCS", "fCS", "CS.chi0", "P0.CS"]


def shape_reads(shape, M):
    """reads[k] = (reads_inputs, [earlier moves read])"""
    R = []
    for k in range(M):
        if shape == "CHW":
            R.append((True, list(range(k))))
        elif shape == "ST":                     # only the last move reads moves: 1 deep move
            R.append((True, list(range(k)) if k == M - 1 else []))
        elif shape == "SR":                     # only move 0 reads inputs: 1 co-deep move with input
            R.append((k == 0, list(range(k))))
        else:
            raise ValueError(shape)
    return R


class Prog:
    def __init__(self, shape, zs, W, menu):
        self.zs, self.W, self.M, self.menu = zs, W, len(zs), menu
        self.reads = shape_reads(shape, self.M)
        # sources: ('x', l) inputs, ('u', k) moves.  Block k (content of move k), block 'R' (readouts).
        self.src = {}
        for k in range(self.M):
            ri, rm = self.reads[k]
            s = ([("x", l) for l in range(W)] if ri else []) + [("u", i) for i in rm]
            self.src[k] = [(s_, w) for s_ in s for w in menu]
        self.srcR = [(s_, w) for s_ in [("x", l) for l in range(W)] + [("u", i) for i in range(self.M)] for w in menu]

    def init(self, rng):
        th = {k: (rng.standard_normal(len(self.src[k])) + 1j * rng.standard_normal(len(self.src[k]))) / np.sqrt(len(self.src[k]))
              for k in range(self.M)}
        th["R"] = np.zeros((self.W, len(self.srcR)), complex)
        return th

    def run(self, G, X, th):
        """X = list of W input probe arrays (N, S).  returns moves u[k] and the term vectors."""
        u = []
        val = lambda s_: X[s_[1]] if s_[0] == "x" else u[s_[1]]
        for k in range(self.M):
            c = sum(t * G.word(w, val(s_)) for t, (s_, w) in zip(th[k], self.src[k]))
            u.append(G.D(self.zs[k], c))
        return u

    def terms_R(self, G, X, u):
        val = lambda s_: X[s_[1]] if s_[0] == "x" else u[s_[1]]
        return [G.word(w, val(s_)) for (s_, w) in self.srcR]

    def outputs(self, G, X, th):
        u = self.run(G, X, th)
        T = self.terms_R(G, X, u)
        return [sum(th["R"][l, p] * T[p] for p in range(len(T))) for l in range(self.W)]


def make_probes(G, W, S, rng):
    if S == 0:      # exact: x_l = (0..I..0) blocks, so the residual is the exact normalised Frobenius error
        X = []
        for l in range(W):
            B = np.zeros((G.N, W * G.N), complex); B[:, l * G.N:(l + 1) * G.N] = np.eye(G.N); X.append(B)
        return X
    X = [(rng.standard_normal((G.N, S)) + 1j * rng.standard_normal((G.N, S))) / np.sqrt(2) for _ in range(W)]
    return X


def target(G, X, l):
    return G.K(X[l])


def rel_res(P, Gs, Xs, th):
    num = den = 0.0
    for G, X in zip(Gs, Xs):
        O = P.outputs(G, X, th)
        for l in range(P.W):
            T = target(G, X, l)
            num += np.linalg.norm(O[l] - T) ** 2; den += np.linalg.norm(T) ** 2
    return np.sqrt(num / den)


def solve_R(P, Gs, Xs, th):
    """readouts: exact LS (each output separately, normalised per f by the target norm)."""
    for l in range(P.W):
        A = []; y = []
        for G, X in zip(Gs, Xs):
            u = P.run(G, X, th)
            T = P.terms_R(G, X, u)
            A.append(np.stack([t.ravel() for t in T], 1)); y.append(target(G, X, l).ravel())
        A = np.concatenate(A); y = np.concatenate(y)
        th["R"][l] = np.linalg.lstsq(A, y, rcond=None)[0]


def solve_block(P, Gs, Xs, th, k):
    """content of move k: outputs are affine-linear in th[k] (with everything else fixed)."""
    n = len(th[k]); cols = []; base = []
    th0 = dict(th); th0[k] = np.zeros(n, complex)
    for G, X in zip(Gs, Xs):
        O0 = P.outputs(G, X, th0)
        base.append(np.concatenate([(target(G, X, l) - O0[l]).ravel() for l in range(P.W)]))
    base = np.concatenate(base)
    O0s = [P.outputs(G, X, th0) for G, X in zip(Gs, Xs)]
    for p in range(n):
        e = np.zeros(n, complex); e[p] = 1
        th1 = dict(th0); th1[k] = e
        col = []
        for (G, X), O0 in zip(zip(Gs, Xs), O0s):
            O1 = P.outputs(G, X, th1)
            col.append(np.concatenate([(O1[l] - O0[l]).ravel() for l in range(P.W)]))
        cols.append(np.concatenate(col))
    A = np.stack(cols, 1)
    lam = 1e-9
    sol = np.linalg.solve(A.conj().T @ A + lam * np.eye(n), A.conj().T @ base)
    th[k] = sol


def main():
    shape = sys.argv[1]; zs = [int(x) for x in sys.argv[2].split(",")]
    FS = [int(x) for x in sys.argv[3].split(",")]; sweeps = int(sys.argv[4]); seed = int(sys.argv[5])
    evalmax = int(sys.argv[6]) if len(sys.argv) > 6 else 8
    W = 2; S = 0
    rng = np.random.default_rng(seed)
    P = Prog(shape, zs, W, MENU)
    Gs = [Game(f) for f in FS]; Xs = [make_probes(G, W, S, rng) for G in Gs]
    th = P.init(rng)
    t0 = time.time()
    solve_R(P, Gs, Xs, th)
    print("# floor.py %s z=%s fit f=%s sweeps=%d seed=%d menu=%d words" % (shape, zs, FS, sweeps, seed, len(MENU)), flush=True)
    for s in range(sweeps):
        for k in range(P.M):
            solve_block(P, Gs, Xs, th, k)
            # gauge: rescale move k's content to unit norm and compensate in every reader of u_k (function unchanged)
            nrm = np.linalg.norm(th[k])
            if nrm > 0:
                th[k] = th[k] / nrm
                for k2 in range(k + 1, P.M):
                    for p, (s_, w) in enumerate(P.src[k2]):
                        if s_ == ("u", k): th[k2][p] *= nrm
                for p, (s_, w) in enumerate(P.srcR):
                    if s_ == ("u", k): th["R"][:, p] *= nrm
            solve_R(P, Gs, Xs, th)
        if s % 5 == 4 or s == sweeps - 1:
            print("  sweep %d fit residual %.5f  max|R| %.2f  (%.0f s)" % (s + 1, rel_res(P, Gs, Xs, th), np.abs(th["R"]).max(), time.time() - t0), flush=True)
    # evaluate the SAME coefficients at f = 2..evalmax on fresh probes
    rng2 = np.random.default_rng(10_000 + seed)
    out = []
    for f in range(2, evalmax + 1):
        G = Game(f); X = make_probes(G, W, 4 if f <= 6 else 2, rng2)
        r = rel_res(P, [G], [X], th)
        out.append((f, r))
        print("  FIXED-a residual at f=%d: %.5f" % (f, r), flush=True)
    np.save("th_%s_%s_%s_%d.npy" % (shape, sys.argv[2], sys.argv[3], seed), np.array([th], dtype=object), allow_pickle=True)
    print("# DONE", flush=True)


if __name__ == "__main__":
    main()
