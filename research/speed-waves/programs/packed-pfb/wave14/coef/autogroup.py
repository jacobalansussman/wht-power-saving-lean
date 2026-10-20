# Predicted quantisation-noise floor of every split of the 7 symmetric tap pairs (t, 13-t) into contiguous runs, with the
# best power-of-two shift per group, for s = 17 and 19 (white-noise model: each raw step adds 4^shift / 12 per tap).
import sys, itertools, numpy as np
N, T = 1024, 14
h = np.load('proto_h_B.npy'); H = h.reshape(T, N); S = h.sum()
pairs = [(p, 13 - p) for p in range(7)]
def lmax(taps): return np.abs(H[list(taps)]).sum(0).max()
res = []
for cuts in itertools.product([0, 1], repeat=6):
    groups, cur = [], [0]
    for i, c in enumerate(cuts):
        if c: groups.append(cur); cur = [i + 1]
        else: cur.append(i + 1)
    groups.append(cur)
    taps = [sum((pairs[p] for p in g), ()) for g in groups]
    Lm = [lmax(t) for t in taps]; n = [len(t) * N for t in taps]
    for s in (17, 19):
        K = ((1 << (s - 1)) - 1) / 128
        best = None
        for g in np.geomspace(100, 1e5, 2000):
            sh = [max(0, int(np.ceil(np.log2(g * l / K - 1e-12)))) for l in Lm]
            noise = np.sqrt(sum(nn * 4.0 ** k for nn, k in zip(n, sh)) / 12) / (g * S)
            if best is None or noise < best[0]: best = (noise, g, sh)
        res.append((s, len(groups), 20 * np.log10(best[0]), best[1], [list(t) for t in taps], best[2]))
for s in (17, 19):
    print(f"s = {s}: best predicted rms noise (dB rel DC) per number of groups")
    for ng in range(1, 8):
        r = min((x for x in res if x[0] == s and x[1] == ng), key=lambda x: x[2])
        print(f"  {ng} groups: {r[2]:6.1f} dB  g={r[3]:7.0f} groups={r[4]} shifts={r[5]}")
