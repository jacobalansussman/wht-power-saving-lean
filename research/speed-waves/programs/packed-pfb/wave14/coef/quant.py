# Coefficient quantisation of the SKA-Low-like prototype (proto_h.npy, from proto.py) for the packed PFB.
# For each case: plain rounding, the best of a scale search, and a greedy integer local search on the worst
# stop-band point (symmetric pairs flipped to their other rounding neighbour; accepted only if the worst
# excess over the float design's mask improves). Constraints:
#   peak:  |c| <= 2^(w-1) - 1                      (one product per DSP field, per-tap extraction)
#   L1s:   every branch k: 127*P_k + 128*N_k <= 2^(s-1) - 1 and 128*P_k + 127*N_k <= 2^(s-1)
#          (P_k, N_k = sums of positive / |negative| coefficients of branch k over its 14 taps):
#          the low field of the packed cascade, sum_t b_t c_t, then fits s bits signed for ANY 8-bit b,
#          so the whole 14-tap cascade stays packed (s = 19: the 27-bit A port with an offset; s = 18 without).
# Prints one line per case and writes coefficient files c_<case>.npy (int64, 14,336 taps).
import sys, json, numpy as np
D = __file__.rsplit('/', 1)[0]
import os
PROTO = os.environ.get('PROTO', 'A')
h = np.load(D + ('/proto_h_w20_400.npy' if PROTO == 'A' else '/proto_h_B.npy')); L = h.size; N = 1024; T = L // N
NFS, NFR = 1 << 19, 1 << 21
fp, fs1, fs2 = 0.5, 0.685, 1.5

def grid(NF):
    f = np.arange(NF // 2 + 1) / NF * N
    return f, f <= fp, (f >= fs1) & (f < fs2), f >= fs2
G = {NF: grid(NF) for NF in (NFS, NFR)}

def resp(c, NF):
    H = np.fft.rfft(c.astype(np.float64), NF)
    return H, 20 * np.log10(np.maximum(np.abs(H) / abs(H[0]), 1e-30))

def metrics(c, NF=NFR):
    f, pb, s1, s2 = G[NF]
    H, R = resp(c, NF)
    rip = (R[pb].max() - R[pb].min()) / 2
    # quantisation error floor: |Hq - g Hfloat| relative to DC, g = least-squares gain
    g = np.dot(c, h) / np.dot(h, h)
    E = np.abs(np.fft.rfft(c - g * h, NF)) / abs(H[0])
    return dict(ripple=rip, stop1=R[s1].max(), stop2=R[s2].max(), med2=float(np.median(R[s2])),
                frac80=float(np.mean(R[s1 | s2] <= -80.0)), frac80_far=float(np.mean(R[s2] <= -80.0)),
                errfloor=20 * np.log10(E.max()))

EDGE = float(os.environ.get('EDGE', '-60.2'))   # requirement 60 dB + 0.2 dB margin at the alias edge
def excess(R, NF):
    # objective: worst stop-band level beyond 1.5 ch (dB), plus 100 x any excess over EDGE in [0.685, 1.5) ch.
    # The point to attack is the region-1 peak while it violates EDGE, else the region-2 peak.
    f, pb, s1, s2 = G[NF]
    v1 = R[s1]; v2 = R[s2]
    i1, i2 = v1.argmax(), v2.argmax()
    pen = max(0.0, v1[i1] - EDGE)
    J = v2[i2] + 100 * pen
    return J, (f[s1][i1] if pen > 0 else f[s2][i2])

GROUP = int(os.environ.get('GROUP', '14'))     # taps per packed cascade group (14 = the whole branch)
def branch_ok(c, s):
    if s is None: return True
    C = c.reshape(T, N)                       # C[t, k] = c[k + N t]
    lim = 1 << (s - 1)
    for g0 in range(0, T, GROUP):
        Cg = C[g0:g0 + GROUP]
        P = np.where(Cg > 0, Cg, 0).sum(0); Nn = np.where(Cg < 0, -Cg, 0).sum(0)
        if not (np.all(127 * P + 128 * Nn <= lim - 1) and np.all(128 * P + 127 * Nn <= lim)): return False
    return True

def scale_for(case_peak, s, g):
    c = np.round(h * g).astype(np.int64)
    if case_peak is not None and np.abs(c).max() > case_peak: return None
    if not branch_ok(c, s): return None
    return c

def gmax(case_peak, s):
    lo, hi = 1.0, float(1 << 20)
    while hi - lo > 1e-3:
        m = (lo + hi) / 2
        if scale_for(case_peak, s, m) is None: hi = m
        else: lo = m
    return lo

def local_search(c, s, peak, iters, x):
    n = np.arange(L // 2); m = L - 1 - n
    H, R = resp(c, NFS); best, fstar = excess(R, NFS)
    for it in range(iters):
        fc = fstar / N                         # cycles/sample
        Hs = np.sum(c * np.exp(-2j * np.pi * fc * np.arange(L)))
        other = c[n] + np.where(x[n] > c[n], 1, -1)       # other rounding neighbour
        ok = np.abs(other - x[n]) < 1
        if peak is not None: ok &= np.abs(other) <= peak
        d = (other - c[n]).astype(np.float64)
        dH = d * (np.exp(-2j * np.pi * fc * n) + np.exp(-2j * np.pi * fc * m))
        gain = np.abs(Hs) - np.abs(Hs + dH)
        gain[~ok] = -np.inf
        order = np.argsort(-gain)[:24]
        improved = False
        for j in order:
            if gain[j] <= 0: break
            c2 = c.copy(); c2[n[j]] = other[j]; c2[m[j]] = other[j]
            if s is not None and not branch_ok(c2, s): continue
            H2, R2 = resp(c2, NFS); e2, f2 = excess(R2, NFS)
            if e2 < best - 1e-4:
                c, best, fstar, improved = c2, e2, f2, True
                break
        if not improved:
            break
    return c, it + 1

cases = [  # name, peak, s, local-search iterations
    ('w18_round', (1 << 17) - 1, None, 0),
    ('w12', 2047, None, 400), ('w11', 1023, None, 400), ('w10', 511, None, 400), ('w9', 255, None, 400),
    ('packed_s19_L1', None, 19, 600),
    ('packed_s18_L1', None, 18, 400),
    ('guarded_w10_g4', 511, None, 0),
]
if __name__ == '__main__':
    only = sys.argv[1:]
    out = []
    print(json.dumps(dict(case='float', **{k: float(v) for k, v in metrics(h * 1e6).items()})))
    print(f"prototype {PROTO}, edge target {EDGE}, group {GROUP}")
    for name, peak, s, iters in cases:
        if only and name not in only: continue
        if peak is not None and s is None and name.endswith('round'):
            g = peak / np.abs(h).max(); c = np.round(h * g).astype(np.int64)
            r = dict(case=name, quantiser='round, peak = full scale', **metrics(c)); print(json.dumps(r)); out.append(r)
            np.save(f"{D}/c_{PROTO}_{name}" + (f"_g{GROUP}" if GROUP != 14 else "") + ".npy", c); continue
        gm = gmax(peak, s)
        # plain rounding at the largest legal scale
        c0 = scale_for(peak, s, gm)
        r = dict(case=name, quantiser='round, largest legal scale', scale=gm, peak=int(np.abs(c0).max()), **metrics(c0))
        print(json.dumps(r)); out.append(r)
        # scale search: 64 scales in [0.80, 1.0] * gm, best worst excess
        best = None
        for g in np.linspace(0.80, 1.0, 64) * gm:
            c = scale_for(peak, s, g)
            if c is None: continue
            e, _ = excess(resp(c, NFS)[1], NFS)
            if best is None or e < best[0]: best = (e, g, c)
        r = dict(case=name, quantiser='scale search (64)', scale=best[1], peak=int(np.abs(best[2]).max()), **metrics(best[2]))
        print(json.dumps(r)); out.append(r)
        if iters:
            x = h * best[1]
            c, its = local_search(best[2], s, peak, iters, x)
            r = dict(case=name, quantiser=f'scale search + greedy local search ({its} steps)', scale=best[1],
                     peak=int(np.abs(c).max()), L1ok_s19=branch_ok(c, 19), **metrics(c))
            print(json.dumps(r)); out.append(r)
        np.save(f"{D}/c_{PROTO}_{name}" + (f"_g{GROUP}" if GROUP != 14 else "") + ".npy", c if iters else best[2])
        sys.stdout.flush()
    with open(f"{D}/quant_{PROTO}_{'_'.join(only) or 'all'}.json", 'w') as fh:
        json.dump(out, fh, indent=1)
