#!/usr/bin/env python3
"""an2-smallnet / frames.py (own code; numpy from our own venv; run with python3 -I).
Raw-array calibration of the frame calculus that every crossing network uses (an-openai sec. 2):
  Phi_U = J diag(i^{q_U(x)}) J,  q_U(x) = Hamming weight of the projection of x onto U (orthonormal basis), mod 4;
  J the normalised Walsh matrix.  Checks, on explicit 2^m x 2^m matrices and random complex arrays:
  (1) Phi_0 = I, Phi_whole = C^{(x)m};  (2) Phi_{U+V} = Phi_U Phi_V for orthogonal U, V;
  (3) Phi_u^2 = R_u (free shift) for an odd line u, so Phi_{u^perp} Phi_u^{-1} = C^{(x)m} R_u;
  (4) a nested step U < V costs dim V - dim U paid moves: Phi_V Phi_U^{-1} = product of one directional kernel per
      orthonormal basis vector of V cap U^perp, each address-conjugate to C on one coordinate;
  (5) the lone pair (W = 2, X on <u>, Y at 0) has a network of 2m paid moves that swaps the data and delivers the full
      transform to both, i.e. paid = W m - 2 sum(dim source) + 2 (descents) with sum = 1, descents = 1."""
import numpy as np, itertools, sys, json, os

def pc(x): return bin(x).count('1')
def dot(x, y): return pc(x & y) & 1

def walsh(m):
    n = 1 << m
    H = np.array([[(-1) ** dot(x, y) for y in range(n)] for x in range(n)], dtype=complex)
    return H / np.sqrt(n)

def phi(basis, m, J=None):
    """basis: list of ints, orthonormal over F_2 (b_i . b_j = delta_ij)."""
    for i, b in enumerate(basis):
        for j, c in enumerate(basis):
            assert dot(b, c) == (1 if i == j else 0), 'basis not orthonormal'
    J = walsh(m) if J is None else J
    n = 1 << m
    d = np.empty(n, dtype=complex)
    for x in range(n):
        p = 0
        for b in basis:
            if dot(x, b): p ^= b
        d[x] = 1j ** (pc(p) % 4)
    return J @ np.diag(d) @ J

def Cm(m):
    C = 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])
    M = np.array([[1]], dtype=complex)
    for _ in range(m): M = np.kron(M, C)
    return M

def shift(u, m):
    n = 1 << m
    R = np.zeros((n, n), dtype=complex)
    for x in range(n): R[x ^ u, x] = 1
    return R

def perp_basis(basis, m):
    """an orthonormal basis of the orthogonal complement of span(basis), when one exists (Gram-Schmidt over F_2
    on the coordinate vectors, keeping odd vectors)."""
    cand = [1 << i for i in range(m)]
    out = []
    span_ok = list(basis)
    # project coordinate vectors off the given orthonormal basis, then orthonormalise greedily
    vecs = []
    for e in cand:
        v = e
        for b in basis:
            if dot(v, b): v ^= b
        vecs.append(v)
    while True:
        changed = False
        for v in vecs:
            w = v
            for b in out:
                if dot(w, b): w ^= b
            if w and dot(w, w) == 1 and all(dot(w, b) == 0 for b in out + list(basis)):
                out.append(w); changed = True
        if not changed: break
    return out

def main():
    rng = np.random.default_rng(15)
    report = {}
    for m in range(1, 6):
        J = walsh(m); n = 1 << m; full = [1 << i for i in range(m)]
        ok1 = np.allclose(phi([], m, J), np.eye(n)) and np.allclose(phi(full, m, J), Cm(m))
        # orthogonal sum: split the coordinates
        ok2 = True
        for k in range(m + 1):
            U = full[:k]; V = full[k:]
            ok2 &= np.allclose(phi(U + V, m, J), phi(U, m, J) @ phi(V, m, J))
        # an odd line u (all odd-weight vectors) and its perp
        ok3 = True; ok4 = True
        for u in range(1, n):
            if pc(u) % 2 == 0: continue
            Pu = phi([u], m, J)
            ok3 &= np.allclose(Pu @ Pu, shift(u, m))
            pb = perp_basis([u], m)
            if len(pb) == m - 1:
                ok3 &= np.allclose(phi(pb, m, J) @ np.linalg.inv(Pu), Cm(m) @ shift(u, m))
                # nested step <u> < whole costs m - 1 directional kernels
                prod = np.eye(n, dtype=complex)
                for z in pb: prod = phi([z], m, J) @ prod
                ok4 &= np.allclose(phi(full, m, J) @ np.linalg.inv(Pu), prod) or np.allclose(Cm(m) @ np.linalg.inv(Pu), prod)
        # (5) the lone pair network on raw arrays, f = 3 columns
        # the port: the heaviest odd vector whose orthogonal complement has an orthonormal basis (needed for a frame)
        u = max((v for v in range(1, n) if pc(v) % 2 == 1 and len(perp_basis([v], m)) == m - 1), key=pc)
        pb = perp_basis([u], m)
        x = rng.normal(size=(n, 3)) + 1j * rng.normal(size=(n, 3))
        y = rng.normal(size=(n, 3)) + 1j * rng.normal(size=(n, 3))
        X, Y = x.copy(), y.copy()          # X on frame <u>, Y on frame 0
        paid = 0
        X = np.linalg.inv(phi([u], m, J)) @ X; paid += 1          # X descends <u> -> 0 (one directional kernel)
        X, Y = Y, X                                                # free gate at the common frame 0: exchange
        X = phi(full, m, J) @ X; paid += m                         # X climbs 0 -> whole: m kernels
        Y = phi(pb, m, J) @ Y; paid += m - 1                       # Y climbs 0 -> u^perp: m - 1 kernels
        okX = np.allclose(X, Cm(m) @ y)
        okY = np.allclose(Y, Cm(m) @ shift(u, m) @ x)
        ok5 = okX and okY and paid == 2 * m and paid == 2 * m - 2 * 1 + 2 * 1
        report[m] = dict(identity_and_whole=bool(ok1), orthogonal_sums=bool(ok2), line_square_is_shift=bool(ok3),
                         nested_step_cost=bool(ok4), lone_pair_2m=bool(ok5), lone_pair_paid=paid, port=u)
        print(f"m={m}: Phi_0=I and Phi_whole=C^m {ok1}; orthogonal sums {ok2}; Phi_u^2=R_u and Phi_perp Phi_u^-1=C^m R_u {ok3};"
              f" nested step = (m-1) kernels {ok4}; lone pair network: outputs exact {okX and okY}, paid {paid} = 2m {ok5}")
        assert ok1 and ok2 and ok3 and ok5
    json.dump(report, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'frames.json'), 'w'), indent=1)

if __name__ == '__main__':
    main()
