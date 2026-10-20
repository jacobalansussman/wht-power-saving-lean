# Stronger optimised quantiser: active-set minimax local search. Starting from the best scale-searched rounding,
# each step evaluates EVERY legal symmetric-pair move (to the other rounding neighbour) on the current worst
# local maxima of the stop-band (64 in [1.5 ch, Nyquist), 16 in [0.685, 1.5) ch) by direct DFT, applies the move
# that lowers the active-set objective most, and keeps it only if the full 2^19-point response confirms it.
#   PROTO=A GROUP=14 quant2.py CASE ITERS      (CASE as in quant.py; the packing bound is re-checked per move)
import os, sys, json, numpy as np
sys.argv, args = sys.argv[:1], sys.argv[1:]
import quant as Q
name, iters = args[0], int(args[1])
peak, s = {c[0]: (c[1], c[2]) for c in Q.cases}[name]
L, N, T = Q.L, Q.N, Q.T
f, pb, s1, s2 = Q.G[Q.NFS]
gm = Q.gmax(peak, s)
best = None
for g in np.linspace(0.80, 1.0, 64) * gm:
    c = Q.scale_for(peak, s, g)
    if c is None: continue
    e, _ = Q.excess(Q.resp(c, Q.NFS)[1], Q.NFS)
    if best is None or e < best[0]: best = (e, g, c)
J, g, c = best; x = Q.h * g
n = np.arange(L // 2); m = L - 1 - n
banned = np.zeros(L // 2, bool)
def active(R):
    out = []
    for mask, k in ((s2, 64), (s1, 16)):
        idx = np.nonzero(mask)[0]; r = R[idx]
        lm = np.nonzero((r[1:-1] >= r[:-2]) & (r[1:-1] >= r[2:]))[0] + 1
        top = lm[np.argsort(-r[lm])[:k]]
        out.append(idx[top])
    return out
acc = 0
for it in range(iters):
    H, R = Q.resp(c, Q.NFS)
    a2, a1 = active(R)
    fa = np.concatenate([a2, a1]); reg1 = np.r_[np.zeros(len(a2), bool), np.ones(len(a1), bool)]
    fr = f[fa] / N
    Hf = H[fa]; H0 = H[0].real
    other = c[n] + np.where(x[n] > c[n], 1, -1)
    ok = (np.abs(other - x[n]) < 1) & ~banned
    if peak is not None: ok &= np.abs(other) <= peak
    d = (other - c[n]).astype(np.float64)
    E = np.exp(-2j * np.pi * np.outer(n, fr)) + np.exp(-2j * np.pi * np.outer(m, fr))   # (pairs, active)
    newH = np.abs(Hf[None, :] + d[:, None] * E) / np.abs(H0 + 2 * d)[:, None]
    dB = 20 * np.log10(np.maximum(newH, 1e-30))
    Jc = dB[:, ~reg1].max(1) + 100 * np.maximum(0, dB[:, reg1].max(1) - Q.EDGE)
    Jc[~ok] = np.inf
    j = int(np.argmin(Jc))
    if not np.isfinite(Jc[j]): break
    c2 = c.copy(); c2[n[j]] = other[j]; c2[m[j]] = other[j]
    if s is not None and not Q.branch_ok(c2, s): banned[j] = True; continue
    e2, _ = Q.excess(Q.resp(c2, Q.NFS)[1], Q.NFS)
    if e2 < J - 1e-4:
        c, J = c2, e2; acc += 1; banned[:] = False
    else:
        banned[j] = True
r = dict(case=name, group=Q.GROUP, proto=Q.PROTO, quantiser=f'scale search + active-set minimax search ({iters} steps, {acc} accepted)',
         scale=g, peak=int(np.abs(c).max()), L1ok_s19=Q.branch_ok(c, 19) if s else None, **{k: float(v) for k, v in Q.metrics(c).items()})
print(json.dumps(r))
tag = f"{Q.PROTO}_{name}" + (f"_g{Q.GROUP}" if Q.GROUP != 14 else "") + os.environ.get("SUFFIX", "")
np.save(f"{Q.D}/c2_{tag}.npy", c)
