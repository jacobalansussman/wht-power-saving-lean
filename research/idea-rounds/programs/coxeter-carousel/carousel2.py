"""Two follow-ups (track wildcard, wht44/ideas, 16 Oct 2026).

(1) Border check: for r = m - 1 carriers, take the least-squares alpha for scalars lam of growing size and run the
    carousel on a stream.  The age residual shrinks like 1/|lam|, but the drain multiplies the read error by lam^m
    every round, so the stream error does not shrink.
(2) Menu of Coxeter-type steps in the block-3 frame game: all c = p o h_1 with p in P = Stab(L_Z) (|P| = 10752).
    Count those whose first m frames are a geodesic L_Z -> L_X (so that c^m lies in P h_full: one paid move per
    round brings every input to the target frame after m rounds), and tabulate the order of c modulo nothing
    (in Sp(6,2)) and the number of distinct frame sequences (Petrie polygons through L_Z).
"""
import itertools
import numpy as np
from carousel import run, vandermonde_solve
from coxeter import sp_h, mul, inv, lag_dist

def border(out, rng):
    for m in (2, 3):
        for scale in (1.0, 10.0, 100.0):
            lam = list(scale * (1 + 0.3 * np.arange(m - 1)))
            alpha, resid, _ = vandermonde_solve(lam, m)
            err, moves, outs = run(m, 1, lam, alpha, T=8, rng=rng)
            out.append(f"m={m} r={m-1} |lam|~{scale:g}: age residual {resid:.1e}; stream max error {err:.1e}")

def P_block3():
    m = 3
    mats = []
    # GL(3,2)
    gl = []
    for bits in itertools.product([0, 1], repeat=9):
        A = np.array(bits, dtype=np.uint8).reshape(3, 3)
        if round(abs(np.linalg.det(A.astype(float)))) % 2 == 1:
            gl.append(A)
    syms = []
    for bits in itertools.product([0, 1], repeat=6):
        S = np.zeros((3, 3), dtype=np.uint8)
        idx = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
        for b, (i, j) in zip(bits, idx):
            S[i, j] = S[j, i] = b
        syms.append(S)
    for A in gl:
        AinvT = inv(A).T.copy()
        for S in syms:
            # block form on (x|z): x' = A x, z' = S A x + A^-T z  (keeps L_Z = {x = 0}; symplectic)
            M = np.zeros((6, 6), dtype=np.uint8)
            M[:3, :3] = A
            M[3:, :3] = mul(S, A)
            M[3:, 3:] = AinvT
            mats.append(M)
    return mats

def is_symplectic(M, m):
    Om = np.zeros((2 * m, 2 * m), dtype=np.int64)
    Om[:m, m:] = np.eye(m, dtype=np.int64)
    Om[m:, :m] = np.eye(m, dtype=np.int64)
    return np.array_equal((M.astype(np.int64).T @ Om @ M.astype(np.int64)) % 2, Om)

def menu(out):
    m = 3
    P = P_block3()
    out.append(f"|P| built = {len(P)}; all symplectic: {all(is_symplectic(M, m) for M in P)}; "
               f"distinct: {len({M.tobytes() for M in P})}")
    h1 = sp_h(m, 0)
    LZ = np.zeros((6, 3), dtype=np.uint8); LZ[3:, :] = np.eye(3, dtype=np.uint8)
    LX = np.zeros((6, 3), dtype=np.uint8); LX[:3, :] = np.eye(3, dtype=np.uint8)
    good, orders, petrie = 0, {}, set()
    good_orders = {}
    for p in P:
        c = mul(p, h1)
        Q = np.eye(6, dtype=np.uint8)
        frames = []
        geo = True
        for j in range(1, m + 1):
            Q = mul(c, Q)
            F = mul(inv(Q), LZ)
            if lag_dist(F, LZ, m) != j:
                geo = False
                break
            frames.append(F)
        # order of c
        R = c.copy(); k = 1
        while not np.array_equal(R, np.eye(6, dtype=np.uint8)):
            R = mul(c, R); k += 1
        orders[k] = orders.get(k, 0) + 1
        if geo and lag_dist(frames[-1], LX, m) == 0:
            good += 1
            good_orders[k] = good_orders.get(k, 0) + 1
            key = tuple(tuple(sorted(tuple(col) for col in np.unique(F.T, axis=0))) for F in frames)
            petrie.add(key)
    out.append(f"block 3: steps c = p h_1 with p in P: {len(P)}; orders of c: {dict(sorted(orders.items()))}")
    out.append(f"block 3: Coxeter-type steps (frames L_Z -> L_X geodesic in 3 rounds): {good}; "
               f"their orders: {dict(sorted(good_orders.items()))}; distinct frame paths (as column sets): {len(petrie)}")

if __name__ == "__main__":
    rng = np.random.default_rng(3)
    out = ["# (1) border check"]
    border(out, rng)
    out.append("# (2) menu of Coxeter-type steps at block 3")
    menu(out)
    txt = "\n".join(out)
    print(txt)
    open("out_carousel2.txt", "w").write(txt + "\n")
