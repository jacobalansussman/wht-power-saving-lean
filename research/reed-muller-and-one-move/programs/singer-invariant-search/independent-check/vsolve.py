#!/usr/bin/env python3
"""vsolve.py (v-rm-singer): independent L-BFGS search of the canonical Singer system (own code; rm-singer's lbf.py,
lma.py, sfield.py NOT imported).  Fibres: one representative per Z_n-orbit of nonzero forms (my vfield.py).
X_Q = Atil Pi_Q (I + Gamma Pi_Q)^-1 Btil, target I_W (x) diag(o_Q); cost = sum over fibre reps of |X_Q - T_Q|_F^2.
Unknowns: complex circulant blocks A (W x M), B (M x W), Gamma (M x M strictly lower), circulant c -> C[a,b]=c[a-b].
Options: --plain (start at the plain schedule; must give 0), --split FILE (W-1 plain + saved 1-array point), noise.
Evidence only.
"""
import argparse, json, time, sys, math
import numpy as np
import torch
sys.path.insert(0, '.')
from vfield import setup, channel_data
import itertools

ap = argparse.ArgumentParser()
ap.add_argument('--m', type=int, required=True)
ap.add_argument('--W', type=int, required=True)
ap.add_argument('--M', type=int, required=True)
ap.add_argument('--coords', default='sd')
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--scale', type=float, default=0.5)
ap.add_argument('--iters', type=int, default=20000)
ap.add_argument('--chunk', type=int, default=200)
ap.add_argument('--maxtime', type=float, default=1200)
ap.add_argument('--plain', action='store_true')
ap.add_argument('--noise', type=float, default=0.0)
ap.add_argument('--save', default='')
ap.add_argument('--init', default='')
ap.add_argument('--tag', default='')
args = ap.parse_args()
torch.manual_seed(args.seed); np.random.seed(args.seed)
torch.set_num_threads(1)
dt = torch.float64

m, W, M = args.m, args.W, args.M
F, basis, C, sig, Sk, pairs = setup(m, args.coords)
n = F.n; N = len(pairs)
# orbit representatives of nonzero forms (rotation classes of the move sequence K; J is covariant)
seen = set(); reps = []
for q in itertools.product([0, 1], repeat=N):
    if not any(q): continue
    K, J = channel_data(F, C, Sk, pairs, q)
    key = min(K[r:] + K[:r] for r in range(n))
    if key in seen: continue
    seen.add(key); reps.append((K, J))
nf = len(reps)
pi = torch.tensor([K for K, J in reps], dtype=dt)          # (nf, n)
oo = torch.tensor([J for K, J in reps], dtype=dt)          # (nf, n)
p = pi.repeat(1, M)                                         # (nf, M n) channel pattern per move
T = torch.zeros(nf, W * n, W * n, dtype=torch.complex128)
for l in range(W):
    T[:, l * n:(l + 1) * n, l * n:(l + 1) * n] = torch.diag_embed(oo).to(torch.complex128)
idx = torch.tensor([[(a - b) % n for b in range(n)] for a in range(n)])
low = torch.tensor([[1.0 if s > t else 0.0 for t in range(M)] for s in range(M)], dtype=dt)


def circ(c):  # c (..., n) complex -> (..., n, n)
    return c[..., idx]


def blocks(X, R, Cc):  # X (R, Cc, n, n) -> (R n, Cc n)
    return X.permute(0, 2, 1, 3).reshape(R * n, Cc * n)


def unpack(theta):
    a = theta[:W * M * n * 2].view(W, M, n, 2)
    b = theta[W * M * n * 2: 2 * W * M * n * 2].view(M, W, n, 2)
    g = theta[2 * W * M * n * 2:].view(M, M, n, 2)
    A = torch.complex(a[..., 0], a[..., 1]); B = torch.complex(b[..., 0], b[..., 1])
    G = torch.complex(g[..., 0], g[..., 1]) * low[:, :, None]
    return A, B, G


# group fibres by |K| and restrict to D = {(t,k): pi_k = 1}: X = Atil[:,D] (I + Gamma[D,D])^-1 Btil[D,:]
groups = []
for wk in sorted(set(int(x) for x in pi.sum(1).tolist())):
    fs = [f for f in range(nf) if int(pi[f].sum()) == wk]
    D = torch.tensor([[t * n + k for t in range(M) for k in range(n) if pi[f, k] > 0] for f in fs])
    groups.append((torch.tensor(fs), D))


def cost(theta, per=False):
    A, B, G = unpack(theta)
    At = blocks(circ(A), W, M); Bt = blocks(circ(B), M, W); Gt = blocks(circ(G), M, M)
    out = []
    for fs, D in groups:
        g = len(fs); d = D.shape[1]
        L = Gt[D[:, :, None], D[:, None, :]] + torch.eye(d, dtype=torch.complex128)[None]
        Bd = Bt[D]                                      # (g, d, Wn)
        Ad = At[:, D].permute(1, 0, 2)                  # (g, Wn, d)
        Y = torch.linalg.solve_triangular(L, Bd, upper=False, unitriangular=True)
        E = Ad @ Y - T[fs]
        out.append((fs, (E.real ** 2 + E.imag ** 2).sum(dim=(1, 2))))
    if per:
        r = torch.zeros(nf, dtype=dt)
        for fs, v in out: r[fs] = v
        return r
    return sum(v.sum() for fs, v in out)


npar = (2 * W * M + M * M) * n * 2


def plain_theta():
    th = torch.zeros(npar, dtype=dt)
    A = torch.zeros(W, M, n, dtype=torch.complex128); B = torch.zeros(M, W, n, dtype=torch.complex128)
    G = torch.zeros(M, M, n, dtype=torch.complex128)
    rs = [F.log[w] for w in basis] if args.coords == 'sd' else list(range(m))
    t = 0
    for l in range(W):
        ts = []
        for r in rs:
            if t >= M: break
            A[l, t, r % n] = 1; B[t, l, (-r) % n] = 1
            for (s0, r0) in ts: pass
            ts.append((t, r)); t += 1
        for i, (ti, ri) in enumerate(ts):
            for (tj, rj) in ts[:i]:
                G[ti, tj, (rj - ri) % n] = 2
    th[:W * M * n * 2] = torch.view_as_real(A).reshape(-1)
    th[W * M * n * 2: 2 * W * M * n * 2] = torch.view_as_real(B).reshape(-1)
    th[2 * W * M * n * 2:] = torch.view_as_real(G).reshape(-1)
    return th


if args.init:
    theta0 = torch.load(args.init)
elif args.plain:
    theta0 = plain_theta()
else:
    theta0 = args.scale * torch.randn(npar, dtype=dt)
if args.noise: theta0 = theta0 + args.noise * torch.randn(npar, dtype=dt)
theta = theta0.clone().requires_grad_(True)
opt = torch.optim.LBFGS([theta], lr=1, max_iter=args.chunk, history_size=50, tolerance_grad=1e-30,
                        tolerance_change=1e-30, line_search_fn='strong_wolfe')
t0 = time.time()
c0 = float(cost(theta.detach()))
hist = [(0, c0)]
print(json.dumps({'start': c0, 'nf': nf, 'n': n, 'npar_complex': npar // 2, 'args': vars(args)}), flush=True)
steps = 0; best = c0; last = c0
while steps < args.iters and time.time() - t0 < args.maxtime:
    def closure():
        opt.zero_grad(); c = cost(theta); c.backward(); return c
    opt.step(closure)
    steps += args.chunk
    c = float(cost(theta.detach()))
    hist.append((steps, c))
    if steps % (args.chunk * 10) == 0:
        print(json.dumps({'step': steps, 'cost': c, 't': round(time.time() - t0, 1),
                          'maxabs': float(theta.detach().abs().max())}), flush=True)
    best = min(best, c)
    if c < 1e-24: break
    if abs(last - c) < 1e-13 * max(1.0, c) and steps > 2000: break
    last = c
per = cost(theta.detach(), per=True)
A, B, G = unpack(theta.detach())
# per-array residual: diagonal blocks
res = {'final': float(per.sum()), 'steps': steps, 't': round(time.time() - t0, 1), 'seed': args.seed,
       'shape': f"{W} in {M}", 'm': m, 'coords': args.coords, 'max_fibre': float(per.max()),
       'maxabs': float(theta.detach().abs().max()), 'tag': args.tag,
       'hist': hist[::max(1, len(hist) // 40)]}
print(json.dumps(res), flush=True)
if args.save: torch.save(theta.detach(), args.save)
