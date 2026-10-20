"""gauge.py (v-t2-pair): gauge-fixed size of a block-4 2-in-7 point.  The system is invariant under v_t -> lam_t v_t
(B_t, Gamma_t. scaled by lam_t; Gamma_.t, A_.t by 1/lam_t).  Minimise the squared parameter norm over the real
log-scalings mu (a convex problem) and report the minimum and the largest entry after balancing, so that coefficient
growth along a gauge orbit is not mistaken for drift.  Usage: gauge.py k1_m2.pt lm_out.pt  (M=7, W=2, d=15)."""
import sys, torch
torch.set_default_dtype(torch.float64)
M, W, d = 7, 2, 15
low = [(s, t) for s in range(M) for t in range(s)]
def parts_from_k1(D):
    x = D['x']; nB = M * W * d
    return D['A'], x[:nB].reshape(M, W, d), x[nB:].reshape(len(low), d)
def parts_from_lm(D):
    pc = D['pc']; nA = W * M * d; nB = M * W * d
    return pc[:nA].reshape(W, M, d), pc[nA:nA + nB].reshape(M, W, d), pc[nA + nB:].reshape(len(low), d)
for f in sys.argv[1:]:
    D = torch.load(f, weights_only=False)
    A, B, G = parts_from_lm(D) if 'pc' in D else parts_from_k1(D)
    a2 = (A.abs() ** 2).sum(dim=(0, 2)); b2 = (B.abs() ** 2).sum(dim=(1, 2)); g2 = (G.abs() ** 2).sum(dim=1)
    mu = torch.zeros(M, requires_grad=True)
    LSs = torch.tensor([s for s, t in low]); LTt = torch.tensor([t for s, t in low])
    opt = torch.optim.LBFGS([mu], max_iter=500, line_search_fn='strong_wolfe')
    def obj():
        return (b2 * torch.exp(2 * mu)).sum() + (a2 * torch.exp(-2 * mu)).sum() + (g2 * torch.exp(2 * (mu[LSs] - mu[LTt]))).sum()
    def clo():
        opt.zero_grad(); o = obj(); o.backward(); return o
    raw = float(obj().detach())
    for _ in range(5): opt.step(clo)
    lam = torch.exp(mu.detach())
    Ab = A / lam[None, :, None]; Bb = B * lam[:, None, None]; Gb = G * (lam[LSs] / lam[LTt])[:, None]
    mx = max(float(Ab.abs().max()), float(Bb.abs().max()), float(Gb.abs().max()))
    print(f"{f}: cost {D.get('cost', D.get('cost_all'))}  raw |p|^2 {raw:.3f}  gauge-fixed |p|^2 {float(obj().detach()):.3f}"
          f"  gauge-fixed max entry {mx:.3f}  raw max entry {max(float(A.abs().max()), float(B.abs().max()), float(G.abs().max())):.3f}")
