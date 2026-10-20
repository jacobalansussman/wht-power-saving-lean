"""scan_frames.py W [maxclasses]: every array uses ANY of its target's 6 frames (plain programs W in 5W). For each
target class (up to k -> 2^s k + c) and each frame choice: which single deletions pass the exact tests (span and
deletion). Fast form: with D(Q) = number of moves where Q = 1 and n(Q) = number of arrays odd at Q, deleting a move
in direction z passes iff D >= n everywhere already and Q(z) = 0 on every tight form (D = n > 0)."""
import sys, itertools, warnings; warnings.filterwarnings("ignore")
import numpy as np
from inst import mons
from packing2 import all_frames, fd, MB
W = int(sys.argv[1]); lim = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9
ML = mons(MB); N = (1 << MB) - 1; P = len(ML)
coeffs = np.array(list(itertools.product([0, 1], repeat=P)), dtype=np.int64)
zb = np.array([[(z >> j) & 1 for j in range(MB)] for z in range(1, N + 1)])
mon = np.stack([zb[:, i] * zb[:, j] for (i, j) in ML], axis=1)
E = ((coeffs @ mon.T) % 2).astype(np.int16)                      # forms x 31 directions (dir z at column z-1)
def parity(k):
    T = fd.target(k); tv = np.array([T[i][j] for (i, j) in ML]); return ((coeffs @ tv) % 2).astype(np.int16)
PAR = {k: parity(k) for k in range(N)}
def canon(ks):
    best = None
    for s in range(MB):
        for c in range(31):
            t = tuple(sorted(((k * 2 ** s) + c) % N for k in ks))
            if best is None or t < best: best = t
    return best
seen = set(); npass_tot = 0; ninst = 0; out = []
for rest in itertools.combinations_with_replacement(range(N), W - 1):
    ks = canon((0,) + rest)
    if ks in seen: continue
    seen.add(ks)
    n = sum(PAR[k] for k in ks)
    frs = [sorted(sorted(F) for F in all_frames(k)) for k in ks]
    nf = [len(f) for f in frs]
    done = set()
    for ch in itertools.product(*[range(n_) for n_ in nf]):
        c = np.zeros(N, dtype=np.int16)
        for l in range(W):
            for z in frs[l][ch[l]]: c[z - 1] += 1
        key = c.tobytes()
        if key in done: continue
        done.add(key); ninst += 1
        D = E @ c
        if (D < n).any(): continue
        tight = (D == n) & (n > 0)
        bad = E[tight].any(axis=0)
        for z in range(1, N + 1):
            if c[z - 1] > 0 and not bad[z - 1]:
                npass_tot += 1
                out.append((ks, ch, z, int(c[z - 1])))
                print(f"PASS ks={ks} frames={ch} delete dir {z} (multiplicity {int(c[z-1])}) tight-before {int(tight.sum())}", flush=True)
    if len(seen) >= lim: break
print(f"W={W}: classes {len(seen)}, distinct plain multisets {ninst}, passing single deletions {npass_tot}", flush=True)
