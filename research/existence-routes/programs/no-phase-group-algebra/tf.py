#!/usr/bin/env python3
"""tf.py (e3-rmbuild): transfer-form solver for the commutative slice WITH a cyclic group of address maps.

Game: A x| Z_n, where A = <dir z, shifts> (abelian, the commutative slice) and Z_n = <P_sigma>, sigma in GL(m,2) of
order n acting on addresses (P_sigma dir(z) P_sigma^-1 = dir(sigma z)).  n = 1 is the commutative slice itself.

Transfer form (PROOF.md of ex2-families, Lemma 1, with ring coefficients): every program with moves t = 1..M in
directions z_t has data Atil in C[Z_n]^{W x M}, Btil in C[Z_n]^{M x W}, G in C[Z_n]^{M x M} strictly lower triangular
(time order), and its input -> delivery map is   T = I + Atil C (I - G C)^{-1} Btil,   C = diag over moves of
(dir(z_t) - 1)  (operator order kept: the ring is not commutative).  The unit condition is T = K I_W in
C[A x| Z_n], i.e. in every induced representation Ind_A chi:  C^n (channels k), P_sigma |k> = |k+1>,
a |k> = chi(P_sigma^-k a P_sigma^k) |k>, so dir(z)|k> = chi(dir(sigma^-k z)) |k>, K|k> = prod_j chi(dir(sigma^-k e_j)).
One character per Z_n-orbit suffices.

Characters used:
  --chars 2 : the order-2 characters chi_Q(dir z) = (-1)^{Q(z)}, Q a binary quadratic form (sector 0, a RELAXATION:
              necessary, not sufficient).
  --chars 4 : all characters chi_q(dir z) = i^{q(z)}, q a Z4-valued quadratic form (the full condition for words whose
              free letters are P_sigma only; shifts are NOT in the coefficient ring here).
Loss = sum over orbit representatives (weighted by orbit size) of || T - K I ||_F^2.
"""
import argparse, itertools, json, math, sys, time
import numpy as np
import torch

PRIM = {2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011}


def mulx(v, m):
    v <<= 1
    if (v >> m) & 1:
        v ^= PRIM[m]
    return v


def singer_table(m, power=1):
    """sigma = multiplication by alpha^power in F_{2^m} (polynomial basis), as a table on 0..2^m-1."""
    N = 1 << m
    t = list(range(N))
    for _ in range(power):
        t = [mulx(v, m) for v in t]
    return t


def order_of(tab):
    N = len(tab); cur = list(range(N)); k = 0
    while True:
        cur = [tab[v] for v in cur]; k += 1
        if cur == list(range(N)):
            return k


def popc(x):
    return bin(x).count("1")


def binary_forms(m):
    pairs = [(j, k) for j in range(m) for k in range(j, m)]
    N = 1 << m
    tabs = []
    for mask in range(1 << len(pairs)):
        row = []
        for z in range(N):
            s = 0
            for b, (j, k) in enumerate(pairs):
                if (mask >> b) & 1 and (z >> j) & 1 and (z >> k) & 1:
                    s ^= 1
            row.append(2 * s)          # value as exponent of i: (-1)^Q = i^(2Q)
        tabs.append(tuple(row))
    return tabs


def z4_forms(m):
    pairs = [(j, k) for j in range(m) for k in range(j + 1, m)]
    N = 1 << m
    tabs = []
    for a in itertools.product(range(4), repeat=m):
        for mask in range(1 << len(pairs)):
            row = []
            for z in range(N):
                s = sum(a[j] for j in range(m) if (z >> j) & 1)
                for b, (j, k) in enumerate(pairs):
                    if (mask >> b) & 1 and (z >> j) & 1 and (z >> k) & 1:
                        s += 2
                row.append(s % 4)
            tabs.append(tuple(row))
    return tabs


def orbit_reps(tabs, sig):
    """orbits of forms under q -> q o sigma; returns list of (rep_table, orbit_size)."""
    index = {t: i for i, t in enumerate(tabs)}
    seen = [False] * len(tabs); reps = []
    N = len(sig)
    for i, t in enumerate(tabs):
        if seen[i]:
            continue
        orb = []; cur = t
        while True:
            j = index[cur]
            if seen[j]:
                break
            seen[j] = True; orb.append(j)
            cur = tuple(cur[sig[z]] for z in range(N))
        reps.append((t, len(orb)))
    return reps


def build_char_data(m, n, sig, dirs, chars):
    """per orbit rep: c (M*n complex: chi(dir) - 1 per (t,k)), tau (n complex: chi(K) per channel)."""
    N = 1 << m
    # sigma^-k tables
    inv = [0] * N
    for v in range(N):
        inv[sig[v]] = v
    sk = [list(range(N))]
    for k in range(1, n):
        sk.append([inv[v] for v in sk[-1]])     # sigma^-k
    tabs = binary_forms(m) if chars == 2 else z4_forms(m)
    reps = orbit_reps(tabs, sig)
    I4 = np.array([1, 1j, -1, -1j])
    C = []; TAU = []; WT = []
    for t, size in reps:
        c = np.zeros(len(dirs) * n, dtype=complex)
        for ti, z in enumerate(dirs):
            for k in range(n):
                c[ti * n + k] = I4[t[sk[k][z]]] - 1
        tau = np.zeros(n, dtype=complex)
        for k in range(n):
            e = sum(t[sk[k][1 << j]] for j in range(m)) % 4
            tau[k] = I4[e]
        C.append(c); TAU.append(tau); WT.append(size)
    return np.array(C), np.array(TAU), np.array(WT, dtype=float)


class Model:
    def __init__(s, W, M, n, scale=0.7, seed=0, dtype=torch.complex128):
        g = torch.Generator().manual_seed(seed)
        def rnd(*shape):
            return (torch.randn(*shape, generator=g, dtype=torch.float64) +
                    1j * torch.randn(*shape, generator=g, dtype=torch.float64)) * scale
        s.W, s.M, s.n = W, M, n
        s.A = rnd(W, M, n).to(dtype).requires_grad_()
        s.B = rnd(M, W, n).to(dtype).requires_grad_()
        s.G = rnd(M, M, n).to(dtype).requires_grad_()
        s.mask = torch.tril(torch.ones(M, M, dtype=torch.float64), -1).to(dtype)[:, :, None]
        a = torch.arange(n)
        s.idx = (a[:, None] - a[None, :]) % n
    def params(s):
        return [s.A, s.B, s.G]
    def big(s, X, r, c):
        # X (r, c, n) -> circulant blocks -> (r*n, c*n)
        Cm = X[:, :, s.idx]                      # (r, c, n, n)
        return Cm.permute(0, 2, 1, 3).reshape(r * s.n, c * s.n)
    def mats(s):
        W, M, n = s.W, s.M, s.n
        return (s.big(s.A, W, M), s.big(s.B, M, W), s.big(s.G * s.mask, M, M))


def loss_fn(model, Cc, TAU, WT, retres=False):
    W, M, n = model.W, model.M, model.n
    Ab, Bb, Gb = model.mats()
    Mn = M * n
    I = torch.eye(Mn, dtype=Ab.dtype)
    Lmat = I[None] - Gb[None] * Cc[:, None, :]          # (Q, Mn, Mn)
    Y = torch.linalg.solve(Lmat, Bb[None].expand(Cc.shape[0], -1, -1))
    X = Ab[None] @ (Cc[:, :, None] * Y)                  # (Q, Wn, Wn) = T - I
    eyeW = torch.eye(W, dtype=Ab.dtype)
    targ = torch.einsum("ab,qk->qakbk", eyeW, TAU - 1).reshape(-1, W * n, W * n) if False else None
    # target: kron(I_W, diag(tau - 1)) with index (l, k)
    D = torch.diag_embed(TAU - 1)                        # (Q, n, n)
    targ = torch.einsum("ab,qij->qaibj", eyeW, D).reshape(-1, W * n, W * n)
    R = X - targ
    per = (R.abs() ** 2).sum(dim=(1, 2))
    L = (per * WT).sum()
    return (L, per) if retres else L


def run(args):
    m, n = args.m, args.n
    sig = singer_table(m, args.power) if n > 1 else list(range(1 << m))
    if n > 1:
        assert order_of(sig) == n, (order_of(sig), n)
    dirs = [int(d) for d in args.dirs.split(",")]
    M = len(dirs)
    Cc, TAU, WT = build_char_data(m, n, sig, dirs, args.chars)
    if args.orbits:
        rng = np.random.default_rng(args.seed + 1000)
        keep = np.union1d(np.sort(rng.choice(len(WT), size=min(args.orbits, len(WT)), replace=False)), [0])
        Cc, TAU, WT = Cc[keep], TAU[keep], WT[keep]
    Cc = torch.tensor(Cc); TAU = torch.tensor(TAU); WT = torch.tensor(WT)
    best = None
    out = []
    for r in range(args.restarts):
        model = Model(args.W, M, n, scale=args.scale, seed=args.seed + r)
        if args.gscale >= 0:
            with torch.no_grad():
                model.G.mul_(args.gscale / args.scale)
        opt = torch.optim.LBFGS(model.params(), lr=1, max_iter=args.iters, history_size=50,
                                line_search_fn="strong_wolfe", tolerance_grad=1e-14, tolerance_change=1e-16)
        def closure():
            opt.zero_grad()
            L = loss_fn(model, Cc, TAU, WT)
            L.backward()
            return L
        t0 = time.time()
        for rnd_ in range(args.rounds):
            opt.step(closure)
            with torch.no_grad():
                L = loss_fn(model, Cc, TAU, WT).item()
            if L < 1e-20:
                break
        rec = dict(restart=r, loss=L, secs=round(time.time() - t0, 1))
        out.append(rec)
        print(json.dumps(rec), flush=True)
        if best is None or L < best[0]:
            best = (L, model)
    L, model = best
    print("BEST", L, flush=True)
    if args.save:
        torch.save(dict(A=model.A.detach(), B=model.B.detach(), G=(model.G * model.mask).detach(), dirs=dirs,
                        m=m, n=n, power=args.power, W=args.W, chars=args.chars, loss=L), args.save)
    return L


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--n", type=int, default=1)
    ap.add_argument("--power", type=int, default=1, help="sigma = alpha^power (Singer power)")
    ap.add_argument("--W", type=int, required=True)
    ap.add_argument("--dirs", required=True, help="comma-separated directions (ints)")
    ap.add_argument("--chars", type=int, default=2)
    ap.add_argument("--restarts", type=int, default=4)
    ap.add_argument("--iters", type=int, default=500)
    ap.add_argument("--rounds", type=int, default=6)
    ap.add_argument("--scale", type=float, default=0.7)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--save", default="")
    ap.add_argument("--orbits", type=int, default=0)
    ap.add_argument("--gscale", type=float, default=-1.0)
    a = ap.parse_args()
    torch.set_num_threads(1)
    run(a)
