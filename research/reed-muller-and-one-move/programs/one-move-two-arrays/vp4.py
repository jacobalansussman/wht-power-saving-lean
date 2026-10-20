#!/usr/bin/env python3
"""vp4.py (t2-pair): vp3.py plus structured starts (--init plain | merge:p | del:p | merge2h:p,q,hpos).
vp3.py (t2-pair): vp2.py plus exact ALS steps on Btil (per-frequency least squares, like the A step) between
L-BFGS chunks (--als k rounds of B-step + A-step per chunk).
vp2.py: vp.py with a fast numpy fibre builder (blocks 5 and 6).
vp.py: VARIABLE-PROJECTION L-BFGS for the canonical Singer no-phase system (rm-singer lbf.py's system and
fibres, unchanged).  The readout Atil enters linearly, and with circulant blocks its least-squares optimum splits by
Fourier frequency of the channel index: for each frequency w and array l an M-dimensional normal-equation solve.
So the objective is  phi(Btil, Gamma) = min_A sum_f ||A Yhat_f - T_f||^2  with  Yhat_f = Pi_f (I + Gamma Pi_f)^-1 Btil,
optimised by L-BFGS over (Btil, Gamma) with the gradient of the envelope (A* held fixed).  Evidence only: a fit is
"found" when the cost on ALL fibres is below --found; it is then saved for an exact step.
Usage: vp.py --m 5 --W 2 --M 9 [--coords sd|poly] [--nfib 120] [--restarts 3] [--seed 1] [--maxtime 900] [--save f.pt]
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
ap.add_argument('--plain', action='store_true')
ap.add_argument('--restarts', type=int, default=3)
ap.add_argument('--chunk', type=int, default=20)
ap.add_argument('--iters', type=int, default=100000)
ap.add_argument('--scale', type=float, default=0.5)
ap.add_argument('--gscale', type=float, default=0.2)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--nfib', type=int, default=0)
ap.add_argument('--ridge', type=float, default=1e-10)
ap.add_argument('--found', type=float, default=1e-18)
ap.add_argument('--stall', type=float, default=1e-6)
ap.add_argument('--maxtime', type=float, default=900.0)
ap.add_argument('--save', default='')
ap.add_argument('--selftest', action='store_true')
ap.add_argument('--als', type=int, default=2)
ap.add_argument('--init', default='rand', help="rand | plain (W chains, M = W m) | merge:p (W=2 chains, move p of both merged) | del:p (chain-2 move p deleted, its successor reads chain 1)")
ap.add_argument('--noise', type=float, default=0.05)
ap.add_argument('--hscale', type=float, default=0.3)
ap.add_argument('--noall', action='store_true')
args = ap.parse_args()
torch.set_default_dtype(torch.float64)
CT = torch.complex128

S = Singer(args.m, args.coords); F = S.F; n = S.n
if args.p == 0:
    d = 1; hexp = [0]
else:
    d = n // math.gcd(n, args.p); hexp = [(args.p * k) % n for k in range(d)]
W, M = args.W, args.M
zexp = [F.log[S.w[t % args.m]] for t in range(M)] if args.plain else [0] * M
z = [F.a_exp(e) for e in zexp]
def build_fibres():
    """fast numpy builder: one representative per orbit of nonzero forms under Q -> Q o (mult by alpha^e), e in hexp."""
    import numpy as np
    m, N = S.m, S.N
    xs = list(range(1, F.q))
    def cmask(c):
        return sum(1 << b for b, (i, j) in enumerate(S.idx) if c[i] and c[j])
    # value of Q at alpha^-k for k in channel order hexp-independent (all n points), as an n-bit integer per form
    pts = [F.a_exp(-k) for k in range(n)]
    pm = np.array([cmask(S.coord[x]) for x in pts], dtype=np.int64)
    Qs = np.arange(1, 1 << N, dtype=np.int64)
    V = np.zeros(len(Qs), dtype=np.uint64)
    for k in range(n):
        bits = np.bitwise_count((Qs & pm[k]).astype(np.uint64)) & 1
        V |= (bits.astype(np.uint64) << np.uint64(k))
    full = np.uint64((1 << n) - 1)
    canon = V.copy()
    for e in hexp[1:]:
        e = np.uint64(e)
        r = ((V >> e) | (V << (np.uint64(n) - e))) & full
        canon = np.minimum(canon, r)
    _, first = np.unique(canon, return_index=True)
    reps = Qs[np.sort(first)]
    tm = []
    for k in range(d):
        T = S.Tk(hexp[k]); tm.append(sum(1 << b for b, (i, j) in enumerate(S.idx) if T[i][j]))
    out = []
    for Qi in reps.tolist():
        v = [(int(V[Qi - 1]) >> k) & 1 for k in range(n)]          # v[k] = Q(alpha^-k)
        Kt = [[v[(hexp[k] - zexp[t]) % n] for k in range(d)] for t in range(M)]
        J = [k for k in range(d) if bin(Qi & tm[k]).count('1') & 1]
        out.append((Kt, J))
    return out
fibres = build_fibres()
Wd, Md = W * d, M * d
PI = torch.tensor([[Kt[t][k] for t in range(M) for k in range(d)] for Kt, J in fibres], dtype=CT)
TT = torch.zeros(len(fibres), Wd, dtype=CT)
for f, (Kt, J) in enumerate(fibres):
    for l in range(W):
        for k in J: TT[f, l * d + k] = 1
print(f"# vp m={args.m} coords={args.coords} d={d} W={W} M={M} plain={args.plain} fibres={len(fibres)}", flush=True)

idx = (torch.arange(d)[:, None] - torch.arange(d)[None, :]) % d
low = [(s, t) for s in range(M) for t in range(s)]
ls = torch.tensor([s for s, t in low], dtype=torch.long); lt = torch.tensor([t for s, t in low], dtype=torch.long)
nB, nG = M * W * d, len(low) * d

def big(X, r, c):
    return X[:, :, idx].permute(0, 2, 1, 3).reshape(r * d, c * d)

def unpackBG(x):
    B = x[:nB].reshape(M, W, d)
    G = torch.zeros(M, M, d, dtype=CT).index_put((ls, lt), x[nB:].reshape(len(low), d))
    return big(B, M, W), big(G, M, M)

# T-hat per fibre: That[f, l, w, c] = sum_i T_f[(l,i), c] e^{-2 pi i w i/d};  T_f = diag(TT[f])
omega = torch.exp(-2j * math.pi * torch.arange(d)[:, None] * torch.arange(d)[None, :] / d).to(CT)   # (w, i)
def That_of(sel):
    T = TT[sel].reshape(-1, W, d)                                  # (g, l, i)
    out = torch.zeros(len(sel), W, d, Wd, dtype=CT)               # (g, l, w, c) with c = (l', i')
    for l in range(W):
        out[:, l, :, l * d:(l + 1) * d] = omega[None, :, :] * T[:, l, None, :]
    return out

GROUPS = {}
def groups_for(sel):
    key = tuple(sel.tolist())
    if key not in GROUPS:
        by = {}
        for j, f in enumerate(key):
            Dc = PI[f].real.nonzero().flatten()
            by.setdefault(len(Dc), []).append((j, f, Dc))
        GROUPS[key] = [(torch.tensor([j for j, _, _ in v]), torch.tensor([f for _, f, _ in v]),
                        torch.stack([D for _, _, D in v]) if k > 0 else None) for k, v in by.items()]
    return GROUPS[key]

THAT = {}
def Yhat_all(x, sel):
    """Yhat-hat: (g, M, w, Wd) complex: FFT over the channel index of Pi (I + Gamma Pi)^-1 B."""
    Bb, Gb = unpackBG(x)
    Y = torch.zeros(len(sel), Md, Wd, dtype=CT)
    for js, fs, Dc in groups_for(sel):
        if Dc is None: continue
        k = Dc.shape[1]
        L = torch.eye(k, dtype=CT)[None] + Gb[Dc[:, :, None], Dc[:, None, :]]
        Yd = torch.linalg.solve(L, Bb[Dc])
        Y = Y.index_put((js[:, None, None].expand(-1, k, Wd), Dc[:, :, None].expand(-1, k, Wd),
                         torch.arange(Wd)[None, None, :].expand(len(js), k, Wd)), Yd)
    Y = Y.reshape(len(sel), M, d, Wd)
    return torch.fft.fft(Y, dim=2)

def solveA(Yh, Th):
    G = torch.einsum('gswc,gtwc->wst', Yh.conj(), Yh)              # (w, M, M)
    h = torch.einsum('gswc,glwc->wls', Yh.conj(), Th)               # (w, W, M)
    tr = torch.diagonal(G, dim1=1, dim2=2).real.mean(dim=1).clamp_min(1e-300)
    G = G + (args.ridge * tr)[:, None, None] * torch.eye(M, dtype=CT)[None]
    Ah = torch.linalg.solve(G[:, None].expand(-1, W, -1, -1), h[..., None])[..., 0]   # (w, W, M)
    return Ah

def cost(x, sel, retA=False):
    key = tuple(sel.tolist())
    if key not in THAT: THAT[key] = That_of(sel)
    Th = THAT[key]
    Yh = Yhat_all(x, sel)
    with torch.no_grad():
        Ah = solveA(Yh.detach(), Th)
    Xh = torch.einsum('wlt,gtwc->glwc', Ah, Yh)
    c = ((Xh - Th).abs() ** 2).sum() / d
    return (c, Ah) if retA else c

def A_time(Ah):   # (w, W, M) -> (W, M, d) circulant first columns
    return torch.fft.ifft(Ah.permute(1, 2, 0), dim=2)

ALL = torch.arange(len(fibres))
if args.selftest:
    g = torch.Generator().manual_seed(3)
    x = torch.complex(torch.randn(nB + nG, generator=g), torch.randn(nB + nG, generator=g)) * 0.3
    sel = ALL[:40]
    c, Ah = cost(x, sel, retA=True)
    A = A_time(Ah); Ab = big(A, W, M); Bb, Gb = unpackBG(x)
    tot = 0.0
    for f in sel.tolist():
        pi = PI[f]; L = torch.eye(Md, dtype=CT) + Gb * pi[None, :]
        X = (Ab * pi[None, :]) @ torch.linalg.solve(L, Bb)
        tot += float(((X - torch.diag(TT[f])).abs() ** 2).sum())
    print(f"# selftest projected {float(c):.10e} direct {tot:.10e} rel {abs(float(c)-tot)/tot:.2e}")
    x2 = Bstep(x, Ah, sel); c2 = float(cost(x2, sel)); print(f"# selftest Bstep cost {float(c):.6e} -> {c2:.6e}")
    # B-step optimality: random perturbation of B must not decrease cost with A fixed
    sys.exit(0)

TB = {}
def TB_of(sel):
    """T-hat for the B step: TBh[g, r, l, w] = sum_c T_f[r, (l, c)] e^{-2 pi i w c/d}, r = (l', i)."""
    T = TT[sel].reshape(-1, W, d)
    out = torch.zeros(len(sel), Wd, W, d, dtype=CT)
    for l in range(W):
        out[:, l * d:(l + 1) * d, l, :] = T[:, l, :, None] * omega.T[None, :, :]
    return out

def Bstep(x, Ah, sel):
    """exact LS for Btil given A (from Ah) and Gamma; returns new x."""
    key = tuple(sel.tolist())
    if key not in TB: TB[key] = TB_of(sel)
    Tb = TB[key]
    A = A_time(Ah); Ab = big(A, W, M); Bb, Gb = unpackBG(x)
    Z = torch.zeros(len(sel), Wd, Md, dtype=CT)
    for js, fs, Dc in groups_for(sel):
        if Dc is None: continue
        k = Dc.shape[1]
        L = torch.eye(k, dtype=CT)[None] + Gb[Dc[:, :, None], Dc[:, None, :]]
        AD = Ab[:, Dc].permute(1, 0, 2)                              # (g, Wd, k)
        Zd = torch.linalg.solve(L.transpose(1, 2), AD.transpose(1, 2)).transpose(1, 2)   # AD L^-1
        Z = Z.index_put((js[:, None, None].expand(-1, Wd, k), torch.arange(Wd)[None, :, None].expand(len(js), Wd, k),
                         Dc[:, None, :].expand(-1, Wd, k)), Zd)
    Zh = torch.fft.fft(Z.reshape(len(sel), Wd, M, d), dim=3)        # (g, r, t, w)
    G = torch.einsum('grsw,grtw->wst', Zh.conj(), Zh)
    h = torch.einsum('grsw,grlw->wls', Zh.conj(), Tb)
    tr = torch.diagonal(G, dim1=1, dim2=2).real.mean(dim=1).clamp_min(1e-300)
    G = G + (args.ridge * tr)[:, None, None] * torch.eye(M, dtype=CT)[None]
    beta = torch.linalg.solve(G[:, None].expand(-1, W, -1, -1), h[..., None])[..., 0]   # (w, l, t) = bhat_tl(-w)
    neg = (-torch.arange(d)) % d
    bhat = beta[neg]                                                # (nu, l, t): bhat_tl(nu)
    b = torch.fft.ifft(bhat.permute(2, 1, 0), dim=2)               # (t, l, d)
    x2 = x.clone(); x2[:nB] = b.reshape(-1)
    return x2

def plain_init(kind, gen):
    """coefficients of W plain chains (directions = the self-dual basis, via Singer powers), then the requested surgery."""
    import cmath
    m = args.m
    a = [(-F.log[S.w[j]]) % n for j in range(m)]          # content factor sigma^{a_j} makes move j act in direction w_j
    chains = []                                             # list of chains: list of move indices
    moves = []                                              # (chain l, position j)
    for l in range(W):
        for j in range(m): moves.append((l, j))
    def circ(e, c=1.0):
        v = torch.zeros(d, dtype=CT); v[e % d] = c; return v
    Mp = len(moves)
    A = torch.zeros(W, Mp, d, dtype=CT); B = torch.zeros(Mp, W, d, dtype=CT); G = torch.zeros(Mp, Mp, d, dtype=CT)
    for t, (l, j) in enumerate(moves):
        B[t, l] = circ(a[j])
        A[l, t] = circ(-a[j])
        for i, (l2, j2) in enumerate(moves[:t]):
            if l2 == l: G[t, i] = circ(a[j] - a[j2], 2.0)
    if kind == 'plain':
        keep = list(range(Mp)); Bn, Gn = B, G
    elif kind.startswith('merge:'):
        p = int(kind.split(':')[1]); t1, t2 = p, m + p           # chain 0 move p and chain 1 move p
        Bn = B.clone(); Gn = G.clone()
        Bn[t1] = B[t1] + B[t2]
        Gn[t1] = G[t1] + G[t2]                                  # merged move reads both chains' earlier moves
        Gn[:, t1] = G[:, t1] + G[:, t2]                         # every later move of either chain reads the merged one
        keep = [t for t in range(Mp) if t != t2]
    elif kind.startswith('merge2h:'):
        # merge chain-0/chain-1 moves at positions p and q, then add one helper move (random content) later
        spec = kind.split(':')[1].split(',')
        if len(spec) == 3: pos = [int(spec[0]), int(spec[1])]; hpos = int(spec[2])
        else: pos = [int(v) for v in spec[0].split('-')]; hpos = int(spec[1])
        Bn = B.clone(); Gn = G.clone(); drop = []
        for pp in pos:
            t1, t2 = pp, m + pp
            Bn[t1] = Bn[t1] + Bn[t2]; Gn[t1] = Gn[t1] + Gn[t2]; Gn[:, t1] = Gn[:, t1] + Gn[:, t2]
            drop.append(t2)
        keep = [t for t in range(Mp) if t not in drop]
        keep_t = torch.tensor(keep)
        Bk = Bn[keep_t]; Gk = Gn[keep_t][:, keep_t]
        order = sorted(range(len(keep)), key=lambda u: (moves[keep[u]][1], moves[keep[u]][0]))
        o = torch.tensor(order); Bk = Bk[o]; Gk = Gk[o][:, o]
        Mk = len(keep); assert Mk + 1 == M, (Mk, M)
        Bf = torch.zeros(M, W, d, dtype=CT); Gf = torch.zeros(M, M, d, dtype=CT)
        old_idx = [i if i < hpos else i + 1 for i in range(Mk)]
        oi = torch.tensor(old_idx)
        Bf[oi] = Bk; Gf[oi[:, None], oi[None, :]] = Gk
        hs = args.hscale
        Bf[hpos] = torch.complex(torch.randn(W, d, generator=gen), torch.randn(W, d, generator=gen)) * hs
        for i in range(hpos):
            Gf[hpos, i] = torch.complex(torch.randn(d, generator=gen), torch.randn(d, generator=gen)) * hs
        for t in range(hpos + 1, M):
            Gf[t, hpos] = torch.complex(torch.randn(d, generator=gen), torch.randn(d, generator=gen)) * hs
        x = torch.cat([Bf.reshape(-1), Gf[ls, lt].reshape(-1)])
        noise = torch.complex(torch.randn(x.shape, generator=gen), torch.randn(x.shape, generator=gen)) * args.noise
        return x + noise
    elif kind.startswith('del:'):
        p = int(kind.split(':')[1]); t2 = m + p
        Bn, Gn = B, G; keep = [t for t in range(Mp) if t != t2]
    else:
        raise SystemExit('bad init')
    keep_t = torch.tensor(keep)
    Bn = Bn[keep_t]; Gn = Gn[keep_t][:, keep_t]
    assert len(keep) == M, (len(keep), M)
    # order: moves are already in a valid time order (chain 0 then chain 1); merged move sits in chain 0's slot
    # lower-triangularity: a merged move t1 < t2 may now read chain-1 moves before t2 that come AFTER t1 in time.
    # Fix the order: put moves in order of (position j, chain l) so that every dependency points backwards.
    order = sorted(range(M), key=lambda u: (moves[keep[u]][1], moves[keep[u]][0]))
    o = torch.tensor(order)
    Bn = Bn[o]; Gn = Gn[o][:, o]
    x = torch.cat([Bn.reshape(-1), Gn[ls, lt].reshape(-1)])
    noise = torch.complex(torch.randn(x.shape, generator=gen), torch.randn(x.shape, generator=gen)) * args.noise
    upper = Gn.clone(); upper[ls, lt] = 0
    if upper.abs().max() > 0: print('# warning: dependencies dropped by ordering', float(upper.abs().max()))
    return x + noise

def run(seed):
    if args.nfib and args.nfib < len(fibres):
        g0 = torch.Generator().manual_seed(1000 + seed)
        sel = torch.randperm(len(fibres), generator=g0)[:args.nfib]
    else:
        sel = ALL
    gen = torch.Generator().manual_seed(seed)
    x0 = torch.complex(torch.randn(nB + nG, generator=gen), torch.randn(nB + nG, generator=gen)) * (args.scale / math.sqrt(2))
    x0[nB:] *= args.gscale / args.scale
    if args.init != 'rand': x0 = plain_init(args.init, gen)
    xr = torch.view_as_real(x0).clone().requires_grad_(True)
    opt = torch.optim.LBFGS([xr], lr=1, max_iter=args.chunk, history_size=100, line_search_fn='strong_wolfe',
                            tolerance_grad=1e-300, tolerance_change=1e-300)
    def closure():
        opt.zero_grad(); c = cost(torch.view_as_complex(opt.param_groups[0]['params'][0]), sel); c.backward(); return c
    t0 = time.time(); prev = float('inf'); it = 0; c = prev
    while it < args.iters and time.time() - t0 < args.maxtime:
        if args.als:
            with torch.no_grad():
                xc = torch.view_as_complex(xr).detach().clone()
                for _ in range(args.als):
                    _, Ah = cost(xc, sel, retA=True)
                    xc = Bstep(xc, Ah, sel)
            xr = torch.view_as_real(xc).clone().requires_grad_(True)
            opt = torch.optim.LBFGS([xr], lr=1, max_iter=args.chunk, history_size=100, line_search_fn='strong_wolfe',
                                    tolerance_grad=1e-300, tolerance_change=1e-300)
        opt.step(closure); it += args.chunk
        xr = opt.param_groups[0]['params'][0]
        with torch.no_grad():
            c = float(cost(torch.view_as_complex(xr), sel))
        print(f"  seed {seed} it {it} cost {c:.6e} t {time.time() - t0:.0f}s", flush=True)
        if not math.isfinite(c) or c < args.found: break
        if prev < float('inf') and prev - c < args.stall * prev: break
        prev = c
    with torch.no_grad():
        x = torch.view_as_complex(xr).detach().clone()
        call, Ah = cost(x, sel if args.noall else ALL, retA=True); call = float(call)
    return x, Ah, c, call, it, time.time() - t0

best = None
for s in range(args.seed, args.seed + args.restarts):
    x, Ah, c, call, its, tt = run(s)
    rec = dict(m=args.m, coords=args.coords, W=W, M=M, plain=args.plain, seed=s, nfib=args.nfib, cost=c, cost_all=call,
               iters=its, secs=round(tt), fibres=len(fibres), found=call < args.found)
    print(json.dumps(rec), flush=True)
    if best is None or call < best[2]: best = (x, Ah, call, s)
if args.save and best is not None:
    torch.save(dict(x=best[0], A=A_time(best[1]), cost_all=best[2], seed=best[3], zexp=zexp, args=vars(args)), args.save)
