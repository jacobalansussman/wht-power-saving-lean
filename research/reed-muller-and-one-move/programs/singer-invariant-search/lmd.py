#!/usr/bin/env python3
"""lmd.py (rm-singer): DENSE Levenberg-Marquardt (Jacobian by torch jacfwd, Marquardt scaling) on the canonical
Singer system.  H = <sigma^p> (order d), moves in field directions z_t, W arrays.  For every H-orbit of nonzero
binary quadratic forms Q (one fibre each), with h_k = sigma^{p k}:
    D = {(t,k): Q(alpha^{-pk} z_t) = 1},   J = {k: Tr(a_Q alpha^{-pk}) = 1}   (linear-part theorem)
    Atil[:, D] (I + Gamma[D, D])^{-1} Btil[D, :] = I_W (x) P_J
with Atil, Btil, Gamma made of d x d circulants (Gamma strictly lower in move order).  This is EXACTLY the no-phase
game with address group H (Theorem A + Lemma C of exist3/e3-rmbuild), so a residual 0 on all fibres is a unit.
d = 1 is the commutative slice (order-2 system of ex2-families).
Prints one JSON line per restart (final residual on all fibres); --save stores the best point (.pt).
"""
import argparse, json, time, sys, math
import torch
from torch.func import jvp, vjp, vmap
sys.path.insert(0, '.')
from sfield import Singer

ap = argparse.ArgumentParser()
ap.add_argument('--m', type=int, required=True)
ap.add_argument('--p', type=int, default=1, help='H = <sigma^p>; p = 0 means H trivial (slice)')
ap.add_argument('--W', type=int, required=True)
ap.add_argument('--M', type=int, required=True)
ap.add_argument('--dirs', default='', help='comma list of exponents e (z_t = alpha^e); default all z_t = 1')
ap.add_argument('--alldirs', action='store_true', help='z_t = alpha^t, t = 0..M-1 (U_m pattern for the slice)')
ap.add_argument('--restarts', type=int, default=2)
ap.add_argument('--iters', type=int, default=200)
ap.add_argument('--cg', type=int, default=40)
ap.add_argument('--chunk', type=int, default=64)
ap.add_argument('--scale', type=float, default=0.5)
ap.add_argument('--gscale', type=float, default=0.3)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--save', default='')
ap.add_argument('--tol', type=float, default=1e-24)
ap.add_argument('--maxtime', type=float, default=3000.0, help='seconds per restart')
ap.add_argument('--coords', default='sd', help="sd: self-dual basis (sigma symmetric); poly: polynomial basis")
ap.add_argument('--plain', action='store_true', help='z_t = self-dual basis vectors, cyclically (plain schedule)')
ap.add_argument('--nfib', type=int, default=0, help='search on a random subset of this many fibres (0 = all)')
args = ap.parse_args()
torch.set_default_dtype(torch.float64)
CT = torch.complex128

S = Singer(args.m, args.coords); F = S.F; n = S.n
if args.p == 0:
    d = 1; hexp = [0]
else:
    d = n // math.gcd(n, args.p); hexp = [(args.p * k) % n for k in range(d)]
W, M = args.W, args.M
if args.dirs:
    zexp = [int(x) for x in args.dirs.split(',')]
elif args.alldirs:
    zexp = list(range(M))
elif args.plain:
    zexp = [F.log[S.w[t % args.m]] for t in range(M)]   # orthonormal basis of the chosen coordinates
else:
    zexp = [0] * M
assert len(zexp) == M
z = [F.a_exp(e) for e in zexp]

# fibres: H-orbits of nonzero forms. Q o h^{-1} : value at x is Q(h^{-1} x); we use values on F to identify forms.
def qvec(Q):
    return tuple(S.Qval(Q, x) for x in range(1, F.q))
seen = set(); fibres = []
for Q in S.forms():
    if not any(Q): continue
    v = qvec(Q)
    if v in seen: continue
    # orbit: Q_h(x) = Q(alpha^{-e} x) for e in hexp
    orb = set()
    for e in hexp:
        orb.add(tuple(S.Qval(Q, F.mul(F.a_exp(-e), x)) for x in range(1, F.q)))
    seen |= orb
    a = S.linpart(Q)
    Dset = [(t, k) for t in range(M) for k in range(d) if S.Qval(Q, F.mul(F.a_exp(-hexp[k]), z[t]))]
    J = [k for k in range(d) if S.pair(Q, S.Tk(hexp[k]))]      # = Tr(a_Q alpha^-pk) in sd coords
    if args.coords == 'sd':
        assert J == [k for k in range(d) if F.tr(F.mul(a, F.a_exp(-hexp[k])))]
    fibres.append((Dset, J, len(orb)))
nodd = sum(1 for f in fibres if f[1])
print(f"# m={args.m} H order d={d} W={W} M={M} dirs={zexp[:20]}{'...' if M > 20 else ''} fibres={len(fibres)} (odd {nodd})",
      flush=True)
# deletion check (exact, necessary): |D| >= W |J|
dele = [(len(D), W * len(J)) for D, J, _ in fibres if len(D) < W * len(J)]
print(f"# deletion violations: {len(dele)} {dele[:5]}", flush=True)

idx = (torch.arange(d)[:, None] - torch.arange(d)[None, :]) % d     # circulant index
mask = torch.tril(torch.ones(M, M), -1)
nA, nB, nG = W * M * d, M * W * d, M * M * d
P = 2 * (nA + nB + nG)

# precompute per-fibre index tensors, padded to one size (batched solve)
def build(fl):
    Dmax = max(1, max(len(D) for D, J, _ in fl))
    cols = torch.zeros(len(fl), Dmax, dtype=torch.long); msk = torch.zeros(len(fl), Dmax)
    T = torch.zeros(len(fl), W * d, W * d, dtype=CT)
    for i, (Dset, J, osz) in enumerate(fl):
        for j, (t, k) in enumerate(Dset):
            cols[i, j] = t * d + k; msk[i, j] = 1
        for l in range(W):
            for k in J:
                T[i, l * d + k, l * d + k] = 1
    return cols, msk.to(CT), T
ALL = build(fibres)
if args.nfib and args.nfib < len(fibres):
    gsel = torch.Generator().manual_seed(1000 + args.seed)
    sel = torch.randperm(len(fibres), generator=gsel)[:args.nfib].tolist()
    SRCH = build([fibres[i] for i in sel])
else:
    SRCH = ALL

def unpack(th):
    re, im = th[:P // 2], th[P // 2:]
    c = torch.complex(re, im)
    A = c[:nA].reshape(W, M, d); B = c[nA:nA + nB].reshape(M, W, d); G = c[nA + nB:].reshape(M, M, d)
    G = G * mask[:, :, None]
    return A, B, G

def big(X, r, cc):
    return X[:, :, idx].permute(0, 2, 1, 3).reshape(r * d, cc * d)

def resid_on(th, dat):
    cols, msk, T = dat
    A, B, G = unpack(th)
    Ab, Bb, Gb = big(A, W, M), big(B, M, W), big(G, M, M)
    AD = Ab[:, cols].permute(1, 0, 2) * msk[:, None, :]
    BD = Bb[cols] * msk[:, :, None]
    GD = Gb[cols[:, :, None], cols[:, None, :]] * (msk[:, :, None] * msk[:, None, :])
    L = torch.eye(cols.shape[1], dtype=CT)[None] + GD
    X = AD @ torch.linalg.solve(L, BD) - T
    r = X.reshape(-1)
    return torch.cat([r.real, r.imag])

def resid(th):
    return resid_on(th, SRCH)

def run(seed):
    g = torch.Generator().manual_seed(seed)
    th = torch.randn(P, generator=g) * args.scale
    gpart = torch.zeros(P, dtype=torch.bool)
    gpart[nA + nB:P // 2] = True; gpart[P // 2 + nA + nB:] = True
    th[gpart] *= args.gscale / args.scale
    lam = 1e-3; t0 = time.time()
    r = resid(th); cost = float(r @ r); hist = [cost]
    def jac(th):
        cols = []
        f = lambda v: jvp(resid, (th,), (v,))[1]
        for c0 in range(0, P, args.chunk):
            c1 = min(P, c0 + args.chunk)
            V = torch.zeros(c1 - c0, P); V[torch.arange(c1 - c0), torch.arange(c0, c1)] = 1
            cols.append(vmap(f)(V).T)
        return torch.cat(cols, dim=1)
    for it in range(args.iters):
        if cost < args.tol or time.time() - t0 > args.maxtime: break
        J = jac(th)
        JTJ = J.T @ J; g_ = J.T @ r; dg = torch.diagonal(JTJ).clone() + 1e-9
        del J
        while True:
            Amat = JTJ + lam * torch.diag(dg)
            try:
                x = torch.linalg.solve(Amat, -g_)
            except Exception:
                lam *= 10; continue
            thn = th + x; rn = resid(thn); cn = float(rn @ rn)
            if cn < cost:
                th, r, cost = thn, rn, cn; lam = max(lam / 5, 1e-15); break
            lam *= 5
            if lam > 1e12: break
        hist.append(cost)
        if lam > 1e12: break
        if it % 5 == 0:
            print(f"  seed {seed} it {it} cost {cost:.6e} lam {lam:.1e} t {time.time() - t0:.0f}s", flush=True)
    return th, cost, hist

best = None
for s in range(args.seed, args.seed + args.restarts):
    th, cost, hist = run(s)
    rall = resid_on(th, ALL); call = float(rall @ rall)
    rec = dict(m=args.m, d=d, W=W, M=M, seed=s, cost=cost, cost_all=call, iters=len(hist) - 1, P=P, fibres=len(fibres),
               searched=int(SRCH[0].shape[0]))
    print(json.dumps(rec), flush=True)
    if best is None or call < best[1]:
        best = (th, call, s)
if args.save and best is not None:
    torch.save(dict(theta=best[0], cost=best[1], seed=best[2], args=vars(args)), args.save)
