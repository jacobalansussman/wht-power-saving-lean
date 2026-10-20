#!/usr/bin/env python3
"""lbfH.py (rm-singer): the L-BFGS search of lbf.py for the no-phase game with an ARBITRARY finite address group H,
given by matrices over F_2 in self-dual coordinates (e3-rmbuild Lemma C, general H):
  channels a in H; P_g|a> = |g a>; a free coefficient c in C[H] acts by R(c)[a, b] = c[a b^-1];
  move t in direction z_t acts on channel a by (-1)^{Q(a^-1 z_t)}; the target's channels are {a : <Q, a^-1 a^-T> = 1}.
Groups: --H singer (the cycle sigma; must reproduce lbf.py), --H gammaL (sigma together with the Frobenius map phi:
x -> x^2, order m n; phi is orthogonal in self-dual coordinates and fixes the direction 1 and the target),
--H frob (phi alone, order m: orthogonal, so Theorem O of w8 applies; a control).
All moves are in direction 1 (allowed when H is transitive on directions, i.e. for singer and gammaL).
--trivial plain: start from the plain schedule on the first m*floor(M/m) moves (exact_plain.py's convention, with
h = sigma^r), for checking the conventions.  Evidence only.
"""
import argparse, json, time, sys, math, itertools
import torch
sys.path.insert(0, '.')
from sfield import Singer

ap = argparse.ArgumentParser()
ap.add_argument('--m', type=int, required=True)
ap.add_argument('--H', default='gammaL')
ap.add_argument('--W', type=int, required=True)
ap.add_argument('--M', type=int, required=True)
ap.add_argument('--restarts', type=int, default=2)
ap.add_argument('--iters', type=int, default=40000)
ap.add_argument('--chunk', type=int, default=200)
ap.add_argument('--scale', type=float, default=0.5)
ap.add_argument('--gscale', type=float, default=0.2)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--found', type=float, default=1e-10)
ap.add_argument('--stall', type=float, default=1e-9)
ap.add_argument('--maxtime', type=float, default=1200.0)
ap.add_argument('--trivial', default='')
ap.add_argument('--noise', type=float, default=0.0)
ap.add_argument('--save', default='')
args = ap.parse_args()
torch.set_default_dtype(torch.float64)
CT = torch.complex128

S = Singer(args.m, 'sd'); F = S.F; n = S.n; m = args.m
mm_, tp = S.mm, S.tp
I = [[int(i == j) for j in range(m)] for i in range(m)]
key = lambda A: tuple(tuple(r) for r in A)
sig = S.sigma
Phi = [[S.coord[F.mul(S.w[j], S.w[j])][i] for j in range(m)] for i in range(m)]   # column j = coords of w_j^2
assert key(mm_(Phi, tp(Phi))) == key(I), "Frobenius not orthogonal"
gens = {'singer': [sig], 'gammaL': [sig, Phi], 'frob': [Phi]}[args.H]
# close the group under multiplication
els = [I]; idxof = {key(I): 0}; frontier = [I]
while frontier:
    nf = []
    for A in frontier:
        for g in gens:
            B = mm_(g, A)
            if key(B) not in idxof:
                idxof[key(B)] = len(els); els.append(B); nf.append(B)
    frontier = nf
d = len(els)
def inv(A):
    for B in els:
        if key(mm_(A, B)) == key(I): return B
    raise ValueError
INV = [idxof[key(inv(A))] for A in els]
MUL = [[idxof[key(mm_(els[a], els[b]))] for b in range(d)] for a in range(d)]
idx = torch.tensor([[MUL[a][INV[b]] for b in range(d)] for a in range(d)])          # a b^-1
def apply(A, x):   # field element x -> A x (coordinates)
    c = S.coord[x]; return S.elem[tuple(sum(A[i][j] & c[j] for j in range(m)) & 1 for i in range(m))]
Tch = [mm_(els[INV[a]], tp(els[INV[a]])) for a in range(d)]                            # a^-1 a^-T
W, M = args.W, args.M
z0 = 1
# fibres: one per H-orbit of nonzero forms (Q ~ Q o g)
seen = set(); fibres = []
nz = list(range(1, F.q))
for Q in S.forms():
    if not any(Q): continue
    v = tuple(S.Qval(Q, x) for x in nz)
    if v in seen: continue
    for A in els:
        seen.add(tuple(S.Qval(Q, apply(A, x)) for x in nz))
    Kc = [S.Qval(Q, apply(els[INV[a]], z0)) for a in range(d)]
    J = [a for a in range(d) if S.pair(Q, Tch[a])]
    fibres.append((Kc, J))
Wd, Md = W * d, M * d
PI = torch.tensor([[Kc[a] for t in range(M) for a in range(d)] for Kc, J in fibres], dtype=CT)
TT = torch.zeros(len(fibres), Wd, dtype=CT)
for f, (Kc, J) in enumerate(fibres):
    for l in range(W):
        for a in J: TT[f, l * d + a] = 1
nodd = int((TT.abs().sum(1) > 0).sum())
print(f"# lbfH m={m} H={args.H} |H|={d} W={W} M={M} fibres={len(fibres)} (target-odd {nodd})", flush=True)

low = [(s, t) for s in range(M) for t in range(s)]
nA, nB, nG = W * M * d, M * W * d, len(low) * d
Pc = nA + nB + nG
ls = torch.tensor([s for s, t in low], dtype=torch.long); lt = torch.tensor([t for s, t in low], dtype=torch.long)

def big(X, r, c):
    return X[:, :, idx].permute(0, 2, 1, 3).reshape(r * d, c * d)

def unpack(x):
    A = x[:nA].reshape(W, M, d); B = x[nA:nA + nB].reshape(M, W, d)
    G = torch.zeros(M, M, d, dtype=CT).index_put((ls, lt), x[nA + nB:].reshape(len(low), d))
    return big(A, W, M), big(B, M, W), big(G, M, M)

GROUPS = {}
def groups():
    if not GROUPS:
        by = {}
        for f in range(len(fibres)):
            Dc = PI[f].real.nonzero().flatten(); by.setdefault(len(Dc), []).append((f, Dc))
        GROUPS['g'] = [(torch.tensor([f for f, _ in v]), torch.stack([D for _, D in v]) if k else None)
                       for k, v in by.items()]
    return GROUPS['g']

def resid_list(x):
    Ab, Bb, Gb = unpack(x); out = []
    for fs, Dc in groups():
        T = torch.diag_embed(TT[fs])
        if Dc is None: out.append((fs, -T)); continue
        k = Dc.shape[1]
        L = torch.eye(k, dtype=CT)[None] + Gb[Dc[:, :, None], Dc[:, None, :]]
        Y = torch.linalg.solve(L, Bb[Dc])
        out.append((fs, Ab[:, Dc].permute(1, 0, 2) @ Y - T))
    return out

def cost(x):
    return sum((R.abs() ** 2).sum() for _, R in resid_list(x))

def rand_x(seed):
    gen = torch.Generator().manual_seed(seed)
    x = torch.complex(torch.randn(Pc, generator=gen), torch.randn(Pc, generator=gen)) * (args.scale / math.sqrt(2))
    x[nA + nB:] *= args.gscale / args.scale
    return x

def plain_x():
    r0 = [F.log[w] for w in S.w]
    hs = lambda r: idxof[key(S.mpow(sig, r))]
    A = torch.zeros(W, M, d, dtype=CT); B = torch.zeros(M, W, d, dtype=CT); G = torch.zeros(M, M, d, dtype=CT)
    for l in range(min(W, M // m)):
        ts = list(range(l * m, (l + 1) * m))
        for i, t in enumerate(ts):
            A[l, t, hs(r0[i])] = 1; B[t, l, hs(-r0[i])] = 1
            for j, s in enumerate(ts):
                if s > t: G[s, t, hs(r0[i] - r0[j])] = 2
    gl = torch.stack([G[s, t] for s, t in low]) if low else torch.zeros(0, d, dtype=CT)
    return torch.cat([A.reshape(-1), B.reshape(-1), gl.reshape(-1)])

def run(seed):
    if args.trivial == 'plain':
        gen = torch.Generator().manual_seed(seed)
        x0 = plain_x() + args.noise * torch.complex(torch.randn(Pc, generator=gen), torch.randn(Pc, generator=gen)) / math.sqrt(2)
        print(f"  seed {seed} plain start cost {float(cost(x0)):.6e} (noise {args.noise})", flush=True)
    else:
        x0 = rand_x(seed)
    xr = torch.view_as_real(x0).clone().requires_grad_(True)
    opt = torch.optim.LBFGS([xr], lr=1, max_iter=args.chunk, history_size=100, line_search_fn='strong_wolfe',
                            tolerance_grad=1e-300, tolerance_change=1e-300)
    def closure():
        opt.zero_grad(); c = cost(torch.view_as_complex(xr)); c.backward(); return c
    t0 = time.time(); prev = float('inf'); it = 0; c = prev
    while it < args.iters and time.time() - t0 < args.maxtime:
        opt.step(closure); it += args.chunk
        with torch.no_grad(): c = float(cost(torch.view_as_complex(xr)))
        print(f"  seed {seed} it {it} cost {c:.6e} t {time.time() - t0:.0f}s", flush=True)
        if not math.isfinite(c) or c < args.found: break
        if prev < float('inf') and prev - c < args.stall * prev: break
        prev = c
    x = torch.view_as_complex(xr).detach().clone()
    return x, c, it, time.time() - t0

best = None
for s in range(args.seed, args.seed + args.restarts):
    x, c, its, tt = run(s)
    with torch.no_grad():
        rows = torch.zeros(W)
        for fs, R in resid_list(x):
            rows += (R.abs() ** 2).sum(2).reshape(len(fs), W, d).sum((0, 2))
    rec = dict(m=m, H=args.H, order=d, W=W, M=M, seed=s, cost=c, iters=its, secs=round(tt), Pc=Pc,
               fibres=len(fibres), found=c < args.found, per_array_rows=[round(float(v), 6) for v in rows])
    print(json.dumps(rec), flush=True)
    if best is None or c < best[1]: best = (x, c, s)
if args.save and best is not None:
    torch.save(dict(x=best[0], cost=best[1], seed=best[2], args=vars(args)), args.save)
