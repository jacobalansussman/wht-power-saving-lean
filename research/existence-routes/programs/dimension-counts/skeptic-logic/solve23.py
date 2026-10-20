# Independent: find an exact 2-in-3 solution in V(x)V at block 2 by damped Gauss-Newton (holomorphic map, complex
# least squares), then the Jacobian rank at the solution vs the generic rank. Own group, own algebra, own Jacobian.
import sys, numpy as np, warnings
warnings.filterwarnings('ignore')
exec(open('indep_core.py').read())
W, s = int(sys.argv[1]), int(sys.argv[2]); seed = int(sys.argv[3])
rng = np.random.default_rng(seed)
slots = [W + j for j in range(s)] + [W + s] * W
P = sum(slots) * dimA
E = []
for l in range(W):
    Xl = np.zeros((d, d * W), complex); Xl[:, l * d:(l + 1) * d] = np.eye(d); E.append(Xl)
tgt = np.concatenate([(np.zeros((d, d * W), complex) + 0) for _ in range(W)], 1)
T = []
for t in range(W):
    X = np.zeros((d, d * W), complex); X[:, t * d:(t + 1) * d] = Kr; T.append(X)
Tv = np.concatenate([x.reshape(-1) for x in T])
def unpack(p):
    out = []; k = 0
    for n_ in slots: out.append(p[k:k + n_ * dimA].reshape(n_, dimA)); k += n_ * dimA
    return out
def forward(p):
    coef = unpack(p); mats = [np.tensordot(c, Abasis, axes=(1, 0)) for c in coef]
    srcs = list(E)
    for j in range(s):
        c = sum(mats[j][i] @ srcs[i] for i in range(len(srcs))); srcs.append(R @ c)
    outs = [sum(mats[s + t][i] @ srcs[i] for i in range(len(srcs))) for t in range(W)]
    return mats, srcs, np.concatenate([o.reshape(-1) for o in outs])
def jac(p):
    mats, srcs, _ = forward(p)
    cols = []
    for kslot in range(len(slots)):
        for i in range(slots[kslot]):
            dsrc = [np.zeros((dimA, d, d * W), complex) for _ in range(W)]
            for j in range(s):
                dc = np.zeros((dimA, d, d * W), complex)
                for ii in range(len(dsrc)): dc += np.matmul(mats[j][ii], dsrc[ii])
                if j == kslot: dc += np.matmul(Abasis, srcs[i])
                dsrc.append(np.matmul(R, dc))
            outs = []
            for t in range(W):
                kk = s + t; do = np.zeros((dimA, d, d * W), complex)
                for ii in range(len(dsrc)): do += np.matmul(mats[kk][ii], dsrc[ii])
                if kk == kslot: do += np.matmul(Abasis, srcs[i])
                outs.append(do.reshape(dimA, -1))
            cols.append(np.concatenate(outs, 1))
    return np.concatenate(cols, 0).T
def rank(J):
    sv = np.linalg.svd(J, compute_uv=False); r = int((sv > 1e-9 * sv[0]).sum())
    return r, sv[r - 1] / sv[0], (sv[r] / sv[0] if r < len(sv) else 0)
best = None
for start in range(int(sys.argv[4]) if len(sys.argv) > 4 else 4):
    p = (rng.normal(size=P) + 1j * rng.normal(size=P)) * 0.3
    lam = 1e-3
    res = forward(p)[2] - Tv; f0 = np.linalg.norm(res)
    for it in range(200):
        J = jac(p)
        # damped GN: (J^H J + lam I) dp = -J^H res, via augmented least squares
        A_ = np.concatenate([J, np.sqrt(lam) * np.eye(P)], 0); b_ = np.concatenate([-res, np.zeros(P)])
        dp = np.linalg.lstsq(A_, b_, rcond=None)[0]
        res2 = forward(p + dp)[2] - Tv
        if np.linalg.norm(res2) < np.linalg.norm(res):
            p = p + dp; res = res2; lam = max(lam / 5, 1e-14)
        else: lam *= 10
        if np.linalg.norm(res) < 1e-12 or lam > 1e8: break
    print(f'start {start}: it {it} residual {np.linalg.norm(res):.3e} max|p| {np.abs(p).max():.3f} lam {lam:.1e}', flush=True)
    if best is None or np.linalg.norm(res) < best[1]: best = (p, np.linalg.norm(res))
    if best[1] < 1e-12: break
p = best[0]
r0 = rank(jac(p)); rg = rank(jac((rng.normal(size=P) + 1j * rng.normal(size=P))))
print('P', P, 'rank at solution', r0, 'generic rank', rg)
