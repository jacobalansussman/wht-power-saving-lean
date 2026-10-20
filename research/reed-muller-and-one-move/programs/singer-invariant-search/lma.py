#!/usr/bin/env python3
"""lma.py (rm-singer): complex Levenberg-Marquardt with an ANALYTIC Jacobian for the canonical Singer system (same
system and fibres as lmc.py / lmd.py).  The residual X_f = Atil_D (I + Gamma_DD)^-1 Btil_D - I_W (x) P_J is holomorphic
in the complex unknowns (A, B, G), so complex Gauss-Newton is exact:  (J^H J + lam diag) dx = -J^H r.
With Y = L^-1 Btil_D and Z = Atil_D L^-1 (L = I + Gamma_DD), dX = dAtil_D Y + Z dBtil_D - Z dGamma_DD Y; every column
of J is a gather of Y, of Z, or a sum of |K| products of them (no solves per column).
--check runs a finite-difference test of the Jacobian and exits.
"""
import argparse, json, time, sys, math
import torch
sys.path.insert(0, '.')
from sfield import Singer

ap = argparse.ArgumentParser()
ap.add_argument('--m', type=int, required=True)
ap.add_argument('--p', type=int, default=1)
ap.add_argument('--coords', default='sd')
ap.add_argument('--W', type=int, required=True)
ap.add_argument('--M', type=int, required=True)
ap.add_argument('--dirs', default='')
ap.add_argument('--alldirs', action='store_true')
ap.add_argument('--plain', action='store_true')
ap.add_argument('--restarts', type=int, default=2)
ap.add_argument('--iters', type=int, default=300)
ap.add_argument('--scale', type=float, default=0.5)
ap.add_argument('--gscale', type=float, default=0.2)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--nfib', type=int, default=0)
ap.add_argument('--save', default='')
ap.add_argument('--tol', type=float, default=1e-26)
ap.add_argument('--maxtime', type=float, default=3000.0)
ap.add_argument('--check', action='store_true')
args = ap.parse_args()
torch.set_default_dtype(torch.float64)
CT = torch.complex128

S = Singer(args.m, args.coords); F = S.F; n = S.n
if args.p == 0:
    d = 1; hexp = [0]
else:
    d = n // math.gcd(n, args.p); hexp = [(args.p * k) % n for k in range(d)]
W, M = args.W, args.M
if args.dirs: zexp = [int(x) for x in args.dirs.split(',')]
elif args.alldirs: zexp = list(range(M))
elif args.plain: zexp = [F.log[S.w[t % args.m]] for t in range(M)]
else: zexp = [0] * M
z = [F.a_exp(e) for e in zexp]

seen = set(); fibres = []
for Q in S.forms():
    if not any(Q): continue
    v = tuple(S.Qval(Q, x) for x in range(1, F.q))
    if v in seen: continue
    seen |= {tuple(S.Qval(Q, F.mul(F.a_exp(-e), x)) for x in range(1, F.q)) for e in hexp}
    Kt = [[S.Qval(Q, F.mul(F.a_exp(-hexp[k]), z[t])) for k in range(d)] for t in range(M)]
    J = [k for k in range(d) if S.pair(Q, S.Tk(hexp[k]))]
    fibres.append((Kt, J))
print(f"# lma m={args.m} coords={args.coords} d={d} W={W} M={M} dirs={zexp[:16]} fibres={len(fibres)}", flush=True)

Wd, Md = W * d, M * d
idx = (torch.arange(d)[:, None] - torch.arange(d)[None, :]) % d
low = [(s, t) for s in range(M) for t in range(s)]          # Gamma[s, t], s > t
nA, nB, nG = W * M * d, M * W * d, len(low) * d
Pc = nA + nB + nG

def prep(fl):
    out = []
    for Kt, J in fl:
        msk = torch.tensor([Kt[t][k] for t in range(M) for k in range(d)], dtype=torch.bool)
        T = torch.zeros(Wd, Wd, dtype=CT)
        for l in range(W):
            for k in J: T[l * d + k, l * d + k] = 1
        out.append((msk, T))
    return out
ALL = prep(fibres)
if args.nfib and args.nfib < len(fibres):
    g0 = torch.Generator().manual_seed(1000 + args.seed)
    SRCH = [ALL[i] for i in torch.randperm(len(fibres), generator=g0)[:args.nfib].tolist()]
else:
    SRCH = ALL

def unpack(x):
    A = x[:nA].reshape(W, M, d); B = x[nA:nA + nB].reshape(M, W, d)
    G = torch.zeros(M, M, d, dtype=CT)
    gl = x[nA + nB:].reshape(len(low), d)
    for i, (s, t) in enumerate(low): G[s, t] = gl[i]
    return A, B, G

def big(X, r, c):
    return X[:, :, idx].permute(0, 2, 1, 3).reshape(r * d, c * d)

def fibre_solve(x, dat):
    A, B, G = unpack(x)
    Ab, Bb, Gb = big(A, W, M), big(B, M, W), big(G, M, M)
    res = []
    for msk, T in dat:
        Dc = msk.nonzero().flatten()
        if len(Dc) == 0:
            res.append((msk, T, None, None, -T)); continue
        L = torch.eye(len(Dc), dtype=CT) + Gb[Dc][:, Dc]
        Y = torch.linalg.solve(L, Bb[Dc])
        Z = torch.linalg.solve(L.T, Ab[:, Dc].T).T
        X = Ab[:, Dc] @ Y - T
        res.append((msk, T, Y, Z, X))
    return res

def cost_of(x, dat):
    return sum(float((X.abs() ** 2).sum()) for *_, X in fibre_solve(x, dat))

kk = torch.arange(d)
def jac_fibre(msk, Y, Z):
    """complex J (Wd*Wd, Pc) for one fibre."""
    Dc = msk.nonzero().flatten()
    Yp = torch.zeros(M, d, Wd, dtype=CT); Zp = torch.zeros(Wd, M, d, dtype=CT)
    if len(Dc):
        Yp.view(Md, Wd)[Dc] = Y; Zp.view(Wd, Md)[:, Dc] = Z
    # J_A[(l1,k1),(l2,k2); l,t,r] = delta(l1,l) Yp[t, k1-r, (l2,k2)]
    sh = (kk[:, None] - kk[None, :]) % d                    # [k1, r] -> k1 - r
    YA = Yp[:, sh, :]                                        # (M, k1, r, Wd)
    JA = torch.zeros(W, d, Wd, W, M, d, dtype=CT)            # (l1,k1),(l2k2) ; l,t,r
    for l in range(W):
        JA[l, :, :, l, :, :] = YA.permute(1, 3, 0, 2)
    JA = JA.reshape(Wd * Wd, nA)
    # J_B[(l1,k1),(l2,k2); t,l,r] = Zp[(l1,k1), t, k2+r] delta(l2,l)
    sh2 = (kk[:, None] + kk[None, :]) % d                    # [k2, r] -> k2 + r
    ZB = Zp[:, :, sh2]                                       # (Wd, M, k2, r)
    JB = torch.zeros(Wd, W, d, M, W, d, dtype=CT)
    for l in range(W):
        JB[:, l, :, :, l, :] = ZB.permute(0, 2, 1, 3)
    JB = JB.reshape(Wd * Wd, nB)
    # J_G[(a),(b); s,t,r] = - sum_{k'} Zp[a, s, k'+r] Yp[t, k', b]
    JG = torch.zeros(Wd, Wd, len(low), d, dtype=CT)
    ss = torch.tensor([s for s, t in low]); tt = torch.tensor([t for s, t in low])
    for r in range(d):
        Zr = Zp[:, :, (kk + r) % d]                          # (Wd, M, k')
        JG[:, :, :, r] = -torch.einsum('aik,ikb->abi', Zr[:, ss, :], Yp[tt])
    JG = JG.reshape(Wd * Wd, nG)
    return torch.cat([JA, JB, JG], dim=1)

def grad_and_JTJ(x, dat):
    sol = fibre_solve(x, dat)
    JTJ = torch.zeros(Pc, Pc, dtype=CT); g = torch.zeros(Pc, dtype=CT); c = 0.0
    for msk, T, Y, Z, X in sol:
        c += float((X.abs() ** 2).sum())
        if Y is None: continue
        Jf = jac_fibre(msk, Y, Z)
        JTJ += Jf.conj().T @ Jf; g += Jf.conj().T @ X.reshape(-1)
    return c, g, JTJ

def rand_x(seed):
    gen = torch.Generator().manual_seed(seed)
    re = torch.randn(Pc, generator=gen); im = torch.randn(Pc, generator=gen)
    x = torch.complex(re, im) * (args.scale / math.sqrt(2))
    x[nA + nB:] *= args.gscale / args.scale
    return x

if args.check:
    x = rand_x(args.seed); v = rand_x(args.seed + 7)
    sol = fibre_solve(x, SRCH[:3]); eps = 1e-6
    sol2 = fibre_solve(x + eps * v, SRCH[:3])
    for (msk, T, Y, Z, X), (*_, X2) in zip(sol, sol2):
        if Y is None: continue
        Jf = jac_fibre(msk, Y, Z)
        fd = ((X2 - X) / eps).reshape(-1); an = Jf @ v
        print("jac check: rel err", float((fd - an).abs().max() / an.abs().max()))
    sys.exit(0)

def run(seed):
    x = rand_x(seed); lam = 1e-3; t0 = time.time()
    cost = cost_of(x, SRCH); it = 0
    for it in range(args.iters):
        if cost < args.tol or time.time() - t0 > args.maxtime: break
        c, g, JTJ = grad_and_JTJ(x, SRCH)
        dg = torch.diagonal(JTJ).real.clone() + 1e-12
        ok = False
        while lam < 1e12:
            dx = torch.linalg.solve(JTJ + lam * torch.diag(dg).to(CT), -g)
            cn = cost_of(x + dx, SRCH)
            if cn < cost:
                x = x + dx; cost = cn; lam = max(lam / 5, 1e-15); ok = True; break
            lam *= 5
        if not ok: break
        if it % 5 == 0:
            print(f"  seed {seed} it {it} cost {cost:.6e} lam {lam:.1e} t {time.time() - t0:.0f}s", flush=True)
    return x, cost, it

best = None
for s in range(args.seed, args.seed + args.restarts):
    x, cost, its = run(s)
    call = cost_of(x, ALL)
    rec = dict(m=args.m, coords=args.coords, d=d, W=W, M=M, seed=s, cost=cost, cost_all=call, iters=its, Pc=Pc,
               fibres=len(fibres), searched=len(SRCH))
    print(json.dumps(rec), flush=True)
    if best is None or call < best[1]: best = (x, call, s)
if args.save and best is not None:
    torch.save(dict(x=best[0], cost_all=best[1], seed=best[2], args=vars(args)), args.save)
