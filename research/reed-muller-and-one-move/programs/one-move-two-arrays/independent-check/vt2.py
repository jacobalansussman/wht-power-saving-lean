#!/usr/bin/env python3
"""vt2.py (v-t2-pair, skeptic): an INDEPENDENT rebuild of the no-phase Singer system of e3-rmbuild Lemma C
(t2-pair's game), written from the description, not from vp.py/sfield.py.

Game: C[Sym_m x| H], H = <sigma> (sigma = multiplication by a primitive alpha of F_{2^m}), moves [z z^T], free part
C[H], target [T] with T = sum_i b_i b_i^T for a basis b of F (T = I in the coordinates dual to b).
For each nonzero quadratic FUNCTION Q on F (one per H-orbit), channels k = 0..n-1:
   Pi[t,k] = Q(alpha^-k z_t),   J[k] = sum_i Q(alpha^-k b_i) mod 2,
   unit  <=>  for every Q:  A Pi (I + G Pi)^-1 B = diag(J) (x) I_W,
A, B, G with n x n circulant blocks (left-regular C[H]), G strictly block-lower-triangular.
Placements: target tensors T up to the normaliser <sigma, Frobenius> (see placements.py); 'sd' = self-dual basis,
'poly' = 1, alpha, ..., alpha^{m-1}; 'orbK' = representative of the K-th orbit (sorted order, see placements list).
Modes: --list (placements), --plain (plain W in Wm check), --eval file.pt (t2-pair's saved point), search (default).
"""
import argparse, itertools, math, sys, time, json
import torch
torch.set_default_dtype(torch.float64)
CT = torch.complex128
ap = argparse.ArgumentParser()
ap.add_argument('--m', type=int, default=4)
ap.add_argument('--W', type=int, default=2)
ap.add_argument('--M', type=int, default=7)
ap.add_argument('--place', default='sd')
ap.add_argument('--list', action='store_true')
ap.add_argument('--plain', action='store_true')
ap.add_argument('--eval', default='')
ap.add_argument('--seeds', default='1')
ap.add_argument('--maxtime', type=float, default=240)
ap.add_argument('--scale', type=float, default=0.5)
ap.add_argument('--hess', action='store_true')
ap.add_argument('--cont', action='store_true')
ap.add_argument('--lm', action='store_true')
ap.add_argument('--init', default='')
ap.add_argument('--noise', type=float, default=0.1)
a = ap.parse_args()
m = a.m; POLY = {3: 0b1011, 4: 0b10011, 5: 0b100101}[m]; q = 1 << m; n = q - 1

def mul(x, y):
    r = 0
    while y:
        if y & 1: r ^= x
        y >>= 1; x <<= 1
        if x >> m: x ^= POLY
    return r
EXP = [1]
for _ in range(n - 1): EXP.append(mul(EXP[-1], 2))
assert len(set(EXP)) == n
LOG = {v: k for k, v in enumerate(EXP)}
def aexp(k): return EXP[k % n]
def tr(x):
    t = 0; y = x
    for _ in range(m): t ^= y; y = mul(y, y)
    return t

# ---- placements: targets T = g g^T (poly coordinates), orbits under N = <sigma, Frobenius>
def bits(x): return [(x >> i) & 1 for i in range(m)]
def matcols(cols):  # matrix with given field elements as columns (poly coords)
    return tuple(tuple(bits(c)[i] for c in cols) for i in range(m))
def mm(A, B): return tuple(tuple(sum(A[i][t] & B[t][j] for t in range(m)) & 1 for j in range(m)) for i in range(m))
def tp(A): return tuple(zip(*A))
def lin(f): return matcols([f(1 << j) for j in range(m)])
SIG = lin(lambda x: mul(x, 2)); FRO = lin(lambda x: mul(x, x))
ID = tuple(tuple(int(i == j) for j in range(m)) for i in range(m))
NG = {ID}; fr = [ID]
while fr:
    new = []
    for g in fr:
        for s in (SIG, FRO):
            h = mm(s, g)
            if h not in NG: NG.add(h); new.append(h)
    fr = new
def rank2(rows):
    rows = [int(''.join(map(str, r)), 2) for r in rows]; r = 0
    basis = []
    for v in rows:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)
def placements():
    rep = {}  # target -> a basis (field elements) b with sum b b^T = T
    for cols in itertools.product(range(1, q), repeat=m):
        if rank2([bits(c) for c in cols]) < m: continue
        G = matcols(cols); T = mm(G, tp(G))
        if T not in rep: rep[T] = cols
    left = set(rep); orbs = []
    while left:
        T0 = min(left); orb = {mm(mm(g, T0), tp(g)) for g in NG}
        orbs.append((T0, orb)); left -= orb
    orbs.sort(key=lambda o: o[0])
    return orbs, rep
def sd_basis():
    for ws in itertools.combinations(range(1, q), m):
        if all(tr(mul(x, y)) == int(x == y) for x in ws for y in ws): return list(ws)

if a.place == 'sd':
    b = sd_basis()
elif a.place == 'poly':
    b = [1 << i for i in range(m)]
else:
    orbs, rep = placements(); b = list(rep[orbs[int(a.place[3:])][0]])
if a.list:
    orbs, rep = placements()
    Tsd = matcols(sd_basis()); Tsd = mm(Tsd, tp(Tsd))
    for K, (T0, orb) in enumerate(orbs):
        tag = ('sd ' if Tsd in orb else '') + ('poly' if ID in orb else '')
        print(f"orb{K}: size {len(orb)} basis {rep[T0]} {tag}")
    sys.exit(0)

# ---- quadratic functions on F: Q(x) = sum_{i<=j} q_ij x_i x_j in poly coordinates (coordinate-free set)
pairs = [(i, j) for i in range(m) for j in range(i, m)]
def Qtab(qv):
    return tuple(sum(c & bits(x)[i] & bits(x)[j] for c, (i, j) in zip(qv, pairs)) & 1 for x in range(q))
W, M = a.W, a.M
if a.plain:
    M = W * m; z = [b[t % m] for t in range(M)]
else:
    z = [1] * M
seen = set(); fib = []
for qv in itertools.product((0, 1), repeat=len(pairs)):
    if not any(qv): continue
    Q = Qtab(qv)
    if Q in seen: continue
    for e in range(n): seen.add(tuple(Q[mul(aexp(-e), x)] for x in range(q)))
    Pi = [[Q[mul(aexp(-k), z[t])] for k in range(n)] for t in range(M)]
    J = [sum(Q[mul(aexp(-k), bi)] for bi in b) & 1 for k in range(n)]
    fib.append((Pi, J))
d = n; F = len(fib); Md, Wd = M * d, W * d
PI = torch.tensor([[P[t][k] for t in range(M) for k in range(d)] for P, J in fib], dtype=CT)       # (F, Md)
TD = torch.tensor([[J[k] for l in range(W) for k in range(d)] for P, J in fib], dtype=CT)         # (F, Wd)
TGT = torch.diag_embed(TD)
print(f"# vt2 m={m} place={a.place} basis={b} W={W} M={M} plain={a.plain} fibres={F} target={float(TD.real.sum()):.1f}",
      flush=True)
cidx = (torch.arange(d)[:, None] - torch.arange(d)[None, :]) % d
def circ(X, r, c):  # X: (r, c, d) first columns -> (r d, c d) block-circulant, C[i,j] = x[(i-j) mod d]
    return X[:, :, cidx].permute(0, 2, 1, 3).reshape(r * d, c * d)
low = [(s, t) for s in range(M) for t in range(s)]
LS = torch.tensor([s for s, t in low]); LT = torch.tensor([t for s, t in low])
def Gfull(g):
    return torch.zeros(M, M, d, dtype=CT).index_put((LS, LT), g)
def Ymat(Bc, Gc):
    Bb = circ(Bc, M, W); Gb = circ(Gfull(Gc), M, M)
    L = torch.eye(Md, dtype=CT)[None] + Gb[None] * PI[:, None, :]
    C = torch.linalg.solve(L, Bb[None].expand(F, -1, -1))
    return PI[:, :, None] * C                                      # (F, Md, Wd)
def direct_cost(Ac, Bc, Gc):
    X = circ(Ac, W, M)[None] @ Ymat(Bc, Gc)
    return ((X - TGT).abs() ** 2).sum()
def proj_A(Y):
    """optimal circulant A by frequency: Xhat[l,w,c] = sum_t ahat[l,t,w] Yhat[t,w,c] (convolution theorem)."""
    Yh = torch.fft.fft(Y.reshape(F, M, d, Wd), dim=2)            # (F, M, w, Wd)
    Th = torch.fft.fft(TGT.reshape(F, W, d, Wd), dim=2)          # (F, W, w, Wd)
    Dm = Yh.permute(2, 0, 3, 1).reshape(d, F * Wd, M)            # (w, rows, M)
    R = Th.permute(2, 0, 3, 1).reshape(d, F * Wd, W)
    Gm = Dm.conj().transpose(1, 2) @ Dm
    Gm = Gm + 1e-13 * Gm.diagonal(dim1=1, dim2=2).real.mean(1).clamp_min(1e-300)[:, None, None] * torch.eye(M, dtype=CT)
    ah = torch.linalg.solve(Gm, Dm.conj().transpose(1, 2) @ R)   # (w, M, W)
    return torch.fft.ifft(ah.permute(2, 1, 0), dim=2)            # (W, M, d) first columns
def proj_cost(Bc, Gc, retA=False):
    Y = Ymat(Bc, Gc)
    with torch.no_grad(): Ac = proj_A(Y.detach())
    X = circ(Ac, W, M)[None] @ Y
    c = ((X - TGT).abs() ** 2).sum()
    return (c, Ac) if retA else c

if a.plain:
    # plain program: array l runs m moves in directions b_0..b_{m-1}; v_t = Pi_t c_t, c_t = x_l - G v (G = +2 on the chain)
    Bc = torch.zeros(M, W, d, dtype=CT); Ac = torch.zeros(W, M, d, dtype=CT); G = torch.zeros(M, M, d, dtype=CT)
    for l in range(W):
        for i in range(m):
            t = l * m + i; Bc[t, l, 0] = 1; Ac[l, t, 0] = 1
            for s in range(l * m, t): G[t, s, 0] = 2
    Gc = G[LS, LT]
    print("plain direct cost", float(direct_cost(Ac, Bc, Gc)), " projected", float(proj_cost(Bc, Gc)))
    sys.exit(0)

if a.eval:
    D = torch.load(a.eval, weights_only=False); x = D['x']; At = D['A']
    nB = M * W * d
    Bc = x[:nB].reshape(M, W, d); Gc = x[nB:].reshape(len(low), d)
    print("saved cost_all", D['cost_all'])
    print("my direct cost with their A", float(direct_cost(At, Bc, Gc)))
    c, Ac = proj_cost(Bc, Gc, retA=True)
    print("my projected cost", float(c), " |A_mine - A_theirs|max", float((Ac - At).abs().max()))
    print("max |B|,|G|,|A|", float(Bc.abs().max()), float(Gc.abs().max()), float(Ac.abs().max()))
    # full gradient wrt all real parameters (A, B, G) at the point
    p = torch.cat([torch.view_as_real(Ac).flatten(), torch.view_as_real(Bc).flatten(), torch.view_as_real(Gc).flatten()])
    nA = Ac.numel() * 2; nBr = Bc.numel() * 2
    def fullc(p):
        A_ = torch.view_as_complex(p[:nA].reshape(W, M, d, 2)); B_ = torch.view_as_complex(p[nA:nA + nBr].reshape(M, W, d, 2))
        G_ = torch.view_as_complex(p[nA + nBr:].reshape(len(low), d, 2))
        return direct_cost(A_, B_, G_)
    pp = p.clone().requires_grad_(True); cc = fullc(pp); g, = torch.autograd.grad(cc, pp)
    print("full cost", float(cc), " |grad|", float(g.norm()), " |grad|/|p|", float(g.norm() / p.norm()))
    print("  |grad_A|", float(g[:nA].norm()), " |grad_B|", float(g[nA:nA + nBr].norm()), " |grad_G|", float(g[nA + nBr:].norm()))
    with torch.no_grad():
        for eps in (1e-3, 1e-4, 1e-5, 1e-6):
            print(f"  step -eps*grad/|grad|, eps={eps:g}: cost change {float(fullc(p - eps * g / g.norm()) - cc):.6e}"
                  f"  (first-order prediction {-eps * float(g.norm()):.6e})")
    if a.lm:
        # Levenberg-Marquardt on the holomorphic residual R(A,B,G) (complex Gauss-Newton), Jacobian by forward mode
        from torch.func import jvp, vmap
        nAc, nBc = Ac.numel(), Bc.numel()
        def resid(pc):
            A_ = pc[:nAc].reshape(W, M, d); B_ = pc[nAc:nAc + nBc].reshape(M, W, d); G_ = pc[nAc + nBc:].reshape(len(low), d)
            return (circ(A_, W, M)[None] @ Ymat(B_, G_) - TGT).reshape(-1)
        pc = torch.cat([Ac.flatten(), Bc.flatten(), Gc.flatten()]).detach().clone(); P = len(pc)
        def jac(pc):
            cols = []
            for c0 in range(0, P, 35):
                E = torch.zeros(min(35, P - c0), P, dtype=CT)
                E[torch.arange(E.shape[0]), torch.arange(c0, c0 + E.shape[0])] = 1
                cols.append(vmap(lambda v: jvp(resid, (pc,), (v,))[1])(E))
            return torch.cat(cols).T                                   # (nres, P)
        # finite-difference check of the Jacobian (holomorphic convention)
        J = jac(pc); v = torch.randn(P, dtype=CT) * 1e-7
        fd = resid(pc + v) - resid(pc); print("  LM: jacobian FD check rel err", float((fd - J @ v).norm() / fd.norm()), flush=True)
        lam = 1e-3; r = resid(pc); c = float((r.abs() ** 2).sum()); t0 = time.time(); k = 0
        while time.time() - t0 < a.maxtime:
            JH = J.conj().T; Hm = JH @ J; gr = JH @ r
            while True:
                dlt = torch.linalg.solve(Hm + lam * torch.diag(Hm.diagonal().real.clamp_min(1e-12)).to(CT), -gr)
                r2 = resid(pc + dlt); c2 = float((r2.abs() ** 2).sum())
                if c2 < c: pc = pc + dlt; r = r2; c = c2; lam = max(lam / 3, 1e-12); break
                lam *= 4
                if lam > 1e12: break
            k += 1
            print(f"  LM step {k} cost {c:.8f} lam {lam:.1e} max|p| {float(pc.abs().max()):.2f} |dlt| {float(dlt.norm()):.2e} t {time.time()-t0:.0f}s", flush=True)
            if lam > 1e12 or c < 1e-20: break
            J = jac(pc)
        torch.save(dict(pc=pc, cost=c), 'lm_out.pt')
    if a.cont:
        # continue from the saved point: plain L-BFGS on ALL parameters (A, B, G), my own loop
        pr = p.clone().requires_grad_(True)
        opt = torch.optim.LBFGS([pr], lr=1, max_iter=50, history_size=100, line_search_fn='strong_wolfe',
                                tolerance_grad=1e-300, tolerance_change=1e-300)
        def clo():
            opt.zero_grad(); c = fullc(pr); c.backward(); return c
        t0 = time.time(); k = 0
        while time.time() - t0 < a.maxtime:
            opt.step(clo); k += 1
            with torch.no_grad(): c = float(fullc(pr))
            gg = torch.autograd.grad(fullc(pr), pr)[0]
            print(f"  cont step {k} cost {c:.8f} |grad| {float(gg.norm()):.3e} max|p| {float(pr.abs().max()):.2f} t {time.time()-t0:.0f}s", flush=True)
    if a.hess:
        t0 = time.time()
        Hm = torch.autograd.functional.hessian(fullc, p)
        ev = torch.linalg.eigvalsh((Hm + Hm.T) / 2)
        print(f"Hessian ({len(p)} params, {time.time()-t0:.0f}s): lowest eigenvalues", [round(float(v), 6) for v in ev[:12]],
              " #neg<-1e-6:", int((ev < -1e-6).sum()), " #|ev|<1e-6:", int((ev.abs() < 1e-6).sum()), " max", float(ev[-1]))
    sys.exit(0)

def canon_start(kind):
    """plain program written in the CANONICAL form (all moves in direction 1): a move in direction b = alpha^e is
    P^e [11^T] P^-e, so B_t = P^-e_t, Gamma_ts = 2 P^(e_s - e_t) on the chain, A_lt = P^e_t (P^j = first column e_j).
    kind 'plainc' (M = W m) or 'merge:p' (W = 2: move p of both arrays merged into one, M = 2m - 1)."""
    e = [LOG[bi] for bi in b]
    if kind == 'plainc':
        moves = [[(l, i)] for l in range(W) for i in range(m)]
    else:
        p_ = int(kind.split(':')[1]); assert W == 2 and M == 2 * m - 1
        moves = [[(0, i)] for i in range(p_)] + [[(1, i)] for i in range(p_)] + [[(0, p_), (1, p_)]] + \
                [[(0, i)] for i in range(p_ + 1, m)] + [[(1, i)] for i in range(p_ + 1, m)]
    assert len(moves) == M
    where = {lt: t for t, mv in enumerate(moves) for lt in mv}
    Bc = torch.zeros(M, W, d, dtype=CT); Ac = torch.zeros(W, M, d, dtype=CT); G = torch.zeros(M, M, d, dtype=CT)
    for t, mv in enumerate(moves):
        for (l, i) in mv:
            Bc[t, l, (-e[i]) % d] += 1; Ac[l, t, e[i] % d] = 1
            for i2 in range(i):
                G[t, where[(l, i2)], (e[i2] - e[i]) % d] += 2
    return Ac, Bc, G[LS, LT]
if a.plain is False and a.seeds.startswith('plainc'):
    Ac, Bc, Gc = canon_start('plainc'); print("canonical plain: direct", float(direct_cost(Ac, Bc, Gc)), "projected", float(proj_cost(Bc, Gc))); sys.exit(0)
def run(seed):
    gen = torch.Generator().manual_seed(seed)
    nB, nG = M * W * d, len(low) * d
    x0 = torch.complex(torch.randn(nB + nG, generator=gen), torch.randn(nB + nG, generator=gen)) * (a.scale / math.sqrt(2))
    x0[nB:] *= 0.4
    if a.init:
        _, B0, G0 = canon_start(a.init)
        x0 = torch.cat([B0.flatten(), G0.flatten()]) + a.noise * x0 / a.scale
    xr = torch.view_as_real(x0).clone().requires_grad_(True)
    opt = torch.optim.LBFGS([xr], lr=1, max_iter=25, history_size=60, line_search_fn='strong_wolfe',
                            tolerance_grad=1e-300, tolerance_change=1e-300)
    def unpack(xc): return xc[:nB].reshape(M, W, d), xc[nB:].reshape(len(low), d)
    def closure():
        opt.zero_grad(); c = proj_cost(*unpack(torch.view_as_complex(xr))); c.backward(); return c
    t0 = time.time(); prev = float('inf'); best = prev; k = 0
    while time.time() - t0 < a.maxtime:
        opt.step(closure); k += 1
        with torch.no_grad(): c = float(proj_cost(*unpack(torch.view_as_complex(xr))))
        if k % 10 == 0: print(f"  seed {seed} step {k} cost {c:.6e} t {time.time()-t0:.0f}s", flush=True)
        if not math.isfinite(c): break
        best = min(best, c)
        if c < 1e-16: break
        if prev - c < 1e-9 * prev: break
        prev = c
    print(json.dumps(dict(m=m, place=a.place, W=W, M=M, seed=seed, cost=c, best=best, secs=round(time.time() - t0),
                          target=float(TD.real.sum()), fibres=F)), flush=True)
for s in [int(v) for v in a.seeds.split(',')]:
    run(s)
