# Graded-group quantiser for the packed PFB (as-4). The 14 taps of every branch are split into groups by TAP INDEX
# (symmetric sets, e.g. the two centre taps, the next four, the outer eight). Each group is one packed DSP cascade
# with its own exact coefficient-sum guard and its own extraction; the group results are added in fabric with a fixed
# LEFT SHIFT per group (power of two, free wiring), so small outer taps get a finer step than the big centre taps.
#   effective coefficient c_eff[k + N t] = c_raw[k + N t] * 2^shift(group of t)
#   legality (per group g, per branch k, on the RAW values): 127*P + 128*N <= 2^(s-1)-1 and 128*P + 127*N <= 2^(s-1)
# SCHEME names: whole (one group = ap-4 'packed'), g7 (ap-4 groups of 7: 0-6, 7-13), and graded schemes below.
# Usage: PROTO=A|B gq.py SCHEME S ITERS [EDGE] [TAG]     -> c3_<PROTO>_<SCHEME>_s<S><TAG>.npy (effective) + json line
import os, sys, json, time, numpy as np
D = os.path.dirname(os.path.abspath(__file__))
N, T = 1024, 14
PROTO = os.environ.get('PROTO', 'B')
h = np.load({'A': f"{D}/proto_h_w20_400.npy", 'B': f"{D}/proto_h_B.npy", 'As': f"{D}/proto_h_As.npy", 'Bs': f"{D}/proto_h_Bs.npy"}[PROTO]); L = h.size
SCHEMES = {
    'whole': [(tuple(range(14)), 0)],
    'g7':    [(tuple(range(7)), 0), (tuple(range(7, 14)), 0)],
    'c2':    [((6, 7), 0), ((0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13), 0)],              # centre pair apart, no shift
    'g3':    [((6, 7), 2), ((4, 5, 8, 9), 1), ((0, 1, 2, 3, 10, 11, 12, 13), 0)],
    'g4':    [((6, 7), 3), ((5, 8), 1), ((3, 4, 9, 10), 1), ((0, 1, 2, 11, 12, 13), 0)],
    'c2s':   [((6, 7), 1), ((0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13), 0)],              # centre pair one step coarser
    'g5':    [((6,), 3), ((7,), 3), ((5, 8), 1), ((3, 4, 9, 10), 1), ((0, 1, 2, 11, 12, 13), 0)],
    'pt':    [((t,), 0) for t in range(14)],                                              # WP486-style per-tap extraction, one scale
    # centre pair UNPACKED (one product per DSP, 2 extra DSPs per pair-lane, 18-bit), the rest graded and packed
    'u4':    [((6, 7), 0, 'U'), ((5, 8), 1), ((3, 4, 9, 10), 1), ((0, 1, 2, 11, 12, 13), 0)],
    'u3':    [((6, 7), 0, 'U'), ((4, 5, 8, 9), 1), ((0, 1, 2, 3, 10, 11, 12, 13), 0)],
    'g4b':   [((6, 7), 3), ((5, 8), 2), ((3, 4, 9, 10), 1), ((0, 1, 2, 11, 12, 13), 0)],   # centre may reach its own bound
    'g4c':   [((6, 7), 3), ((5, 8), 2), ((3, 4, 9, 10), 2), ((0, 1, 2, 11, 12, 13), 1)],
    'g3b':   [((6, 7), 2), ((5, 8), 0), ((0, 1, 2, 3, 4, 9, 10, 11, 12, 13), 0)],
}
fp, fs1, fs2 = 0.5, 0.685, 1.5
NFS, NFR = 1 << 19, 1 << 21
def grid(NF):
    f = np.arange(NF // 2 + 1) / NF * N
    return f, f <= fp, (f >= fs1) & (f < fs2), f >= fs2
G = {NF: grid(NF) for NF in (NFS, NFR)}
def resp(c, NF):
    H = np.fft.rfft(c.astype(np.float64), NF)
    return H, 20 * np.log10(np.maximum(np.abs(H) / abs(H[0]), 1e-30))
def metrics(c, NF=NFR):
    f, pb, s1, s2 = G[NF]; H, R = resp(c, NF)
    return dict(ripple=float((R[pb].max() - R[pb].min()) / 2), stop1=float(R[s1].max()), stop2=float(R[s2].max()),
                med2=float(np.median(R[s2])), frac80=float(np.mean(R[s1 | s2] <= -80.0)), frac80_far=float(np.mean(R[s2] <= -80.0)),
                f_stop2=float(f[s2][R[s2].argmax()]))

UNP = set()
def setup(scheme):
    sh = np.zeros(T, int); gid = np.zeros(T, int)
    global UNP; UNP = set()
    for j, e in enumerate(SCHEMES[scheme]):
        taps, s_ = e[0], e[1]
        if len(e) > 2: UNP.add(j)
        for t in taps: sh[t] = s_; gid[t] = j
    return sh, gid
def raw_of(ceff, sh):   # ceff (L,) -> raw (T, N)
    return ceff.reshape(T, N) >> sh[:, None]
def legal(raw, gid, s):
    lim = 1 << (s - 1)
    for j in range(gid.max() + 1):
        if j in UNP:
            if np.abs(raw[gid == j]).max() > (1 << 17) - 1: return False
            continue
        R = raw[gid == j]; P = np.where(R > 0, R, 0).sum(0); Nn = np.where(R < 0, -R, 0).sum(0)
        if not (np.all(127 * P + 128 * Nn <= lim - 1) and np.all(128 * P + 127 * Nn <= lim)): return False
    return True
def quant(g, sh):       # effective coefficients at effective scale g
    step = (1 << sh)[:, None].astype(float)
    raw = np.round(h.reshape(T, N) * g / step).astype(np.int64)
    return (raw << sh[:, None]).reshape(-1), raw
def gmax(sh, gid, s):
    lo, hi = 1.0, float(1 << 24)
    while hi / lo > 1 + 1e-6:
        m = (lo * hi) ** 0.5
        if legal(quant(m, sh)[1], gid, s): lo = m
        else: hi = m
    return lo
def excess(R, NF, EDGE):
    f, pb, s1, s2 = G[NF]; v1 = R[s1]; v2 = R[s2]
    rip = (R[pb].max() - R[pb].min()) / 2
    return v2.max() + 100 * max(0.0, v1.max() - EDGE) + 100 * max(0.0, rip - 0.195)

if __name__ == '__main__':
    scheme, s, iters = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    EDGE = float(sys.argv[4]) if len(sys.argv) > 4 else -60.5
    TAG = sys.argv[5] if len(sys.argv) > 5 else ''
    sh, gid = setup(scheme); t0 = time.time()
    gm = gmax(sh, gid, s)
    best = None
    for g in np.linspace(0.80, 1.0, 48) * gm:
        c, raw = quant(g, sh)
        if not legal(raw, gid, s): continue
        e = excess(resp(c, NFS)[1], NFS, EDGE)
        if best is None or e < best[0]: best = (e, g, c)
    J, g, c = best
    r0 = dict(stage='scale search', **metrics(c))
    print(json.dumps(dict(proto=PROTO, scheme=scheme, s=s, gmax=gm, g=g, **r0)), flush=True)
    # dithered rounding: NDITH tries of round(x/step + d), d uniform in [-DLEV, DLEV] (symmetric pairs), keep the best
    ND = int(os.environ.get('NDITH', '0')); DL = float(os.environ.get('DLEV', '0.3'))
    rng = np.random.default_rng(int(os.environ.get('SEED', '1')))
    step = (1 << sh)[:, None].astype(float); X = h.reshape(T, N) * g / step
    for i in range(ND):
        d = rng.uniform(-DL, DL, L // 2); d = np.r_[d, d[::-1]].reshape(T, N)
        raw = np.round(X + d).astype(np.int64)
        if not legal(raw, gid, s): continue
        c2 = (raw << sh[:, None]).reshape(-1)
        e2 = excess(resp(c2, NFS)[1], NFS, EDGE)
        if e2 < J: J, c = e2, c2
    if ND: print(json.dumps(dict(stage=f'dither {ND} x {DL}', J=J, **metrics(c))), flush=True)
    if os.environ.get('INIT'):   # continue from a saved coefficient file (a checkpoint) at scale G0
        g = float(os.environ['G0']); c = np.load(os.environ['INIT']); J = excess(resp(c, NFS)[1], NFS, EDGE)
        print(json.dumps(dict(stage='init from ' + os.environ['INIT'], J=J)), flush=True)
    # active-set minimax on symmetric pairs, moves of +-1 raw step (= +-2^shift effective) toward the other neighbour
    f, pb, s1, s2 = G[NFS]
    n = np.arange(L // 2); m = L - 1 - n; tn = n // N; stepn = (1 << sh[tn]).astype(np.int64)
    x = h * g
    banned = np.zeros(L // 2, bool); acc = 0
    def active(R):
        out = []
        for mask, k in ((s2, int(os.environ.get('KA', '64'))), (s1, 16)):
            idx = np.nonzero(mask)[0]; r = R[idx]
            lm = np.nonzero((r[1:-1] >= r[:-2]) & (r[1:-1] >= r[2:]))[0] + 1
            out.append(idx[lm[np.argsort(-r[lm])[:k]]])
        return out
    ck = os.environ.get('CKPT', '')
    for it in range(iters):
        H, R = resp(c, NFS)
        a2, a1 = active(R); fa = np.concatenate([a2, a1]); reg1 = np.r_[np.zeros(len(a2), bool), np.ones(len(a1), bool)]
        fr = f[fa] / N; Hf = H[fa]; H0 = H[0].real
        other = c[n] + np.where(x[n] > c[n], 1, -1) * stepn
        if os.environ.get('WIDE') and it % 2 == 1:   # odd steps: the move in the OTHER direction (|error| up to 1.5 steps)
            other = c[n] - np.where(x[n] > c[n], 1, -1) * stepn
            ok = (np.abs(other - x[n]) < 1.5 * stepn) & ~banned
        else:
            ok = (np.abs(other - x[n]) < stepn) & ~banned
        d = (other - c[n]).astype(np.float64)
        E = np.exp(-2j * np.pi * np.outer(n, fr)) + np.exp(-2j * np.pi * np.outer(m, fr))
        dB = 20 * np.log10(np.maximum(np.abs(Hf[None, :] + d[:, None] * E) / np.abs(H0 + 2 * d)[:, None], 1e-30))
        Jc = dB[:, ~reg1].max(1) + 100 * np.maximum(0, dB[:, reg1].max(1) - EDGE)
        Jc[~ok] = np.inf
        j = int(np.argmin(Jc))
        if not np.isfinite(Jc[j]): break
        MULTI = int(os.environ.get('MULTI', '1'))
        if MULTI > 1:   # try the MULTI best predicted moves together first (one full check), then the best alone
            top = np.argsort(Jc)[:MULTI]; top = top[np.isfinite(Jc[top])]
            c2 = c.copy(); c2[n[top]] = other[top]; c2[m[top]] = other[top]
            if legal(raw_of(c2, sh), gid, s):
                e2 = excess(resp(c2, NFS)[1], NFS, EDGE)
                if e2 < J - 1e-4:
                    c, J = c2, e2; acc += len(top); banned[:] = False
                    if ck and (it + 1) % 50 == 0:
                        np.save(f"{D}/ckpt_{PROTO}_{scheme}_s{s}{TAG}.npy", c)
                        print(json.dumps(dict(ckpt=it + 1, acc=acc, J=J, t=round(time.time() - t0))), flush=True)
                    continue
        c2 = c.copy(); c2[n[j]] = other[j]; c2[m[j]] = other[j]
        if not legal(raw_of(c2, sh), gid, s): banned[j] = True; continue
        e2 = excess(resp(c2, NFS)[1], NFS, EDGE)
        if e2 < J - 1e-4: c, J = c2, e2; acc += 1; banned[:] = False
        else: banned[j] = True
        if ck and (it + 1) % 50 == 0:
            np.save(f"{D}/ckpt_{PROTO}_{scheme}_s{s}{TAG}.npy", c)
            print(json.dumps(dict(ckpt=it + 1, acc=acc, J=J, t=round(time.time() - t0))), flush=True)
    raw = raw_of(c, sh)
    r = dict(proto=PROTO, scheme=scheme, s=s, edge_target=EDGE, g=g, iters=iters, accepted=acc, legal=legal(raw, gid, s),
             raw_peak=[int(np.abs(raw[gid == j]).max()) for j in range(gid.max() + 1)], shifts=[int(x_) for x_ in sh],
             seconds=round(time.time() - t0), **metrics(c))
    print(json.dumps(r), flush=True)
    np.save(f"{D}/c3_{PROTO}_{scheme}_s{s}{TAG}.npy", c)
