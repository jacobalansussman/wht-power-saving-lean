"""scan3.py m [W=3]: my own enumeration of 'plain minus one move' at block m (rm2-solve section 5, from the
definitions).

Target multisets {k_1..k_W} up to k -> 2^s k + c (mod n); each array on any frame of its target (W frames, a
multiset of mW directions); distinct multisets per class; then every single deletion of one copy of a direction,
kept when span and deletion hold (|D_Q| >= r_Q for every form Q).  Prints one line per passing deletion and the
totals.  Fast form (as the definitions give): after deletion of direction z, |D'_Q| = |D_Q| - Q(z), so a deletion
passes iff the plain multiset passes and no tight form (|D_Q| = r_Q >= 1) has Q(z) = 1.
"""
import sys, itertools
import numpy as np
from f2q import Field, all_frames, check_field, popc

m = int(sys.argv[1]); W = int(sys.argv[2]) if len(sys.argv) > 2 else 3
F = Field(m); check_field(F); n = F.n; P = F.P; NQ = 1 << P
Qs = np.arange(NQ, dtype=np.int64)
def evalmask(mask):
    return np.array([popc(int(q) & mask) & 1 for q in range(NQ)], dtype=np.int16)
EV = np.stack([evalmask(F.dir_mask(z)) for z in range(1, n + 1)], axis=1)   # NQ x n, column z-1
PAR = {k: evalmask(F.sym_mask(F.target(k))) for k in range(n)}
FR = {k: all_frames(F, k) for k in range(n)}
print(f"block {m}: self-dual basis {F.basis}; frames per target {sorted(set(len(v) for v in FR.values()))}", flush=True)

def canon(ks):
    best = None
    for s in range(m):
        for c in range(n):
            t = tuple(sorted(((k << s) + c) % n for k in ks))
            if best is None or t < best: best = t
    return best

classes = sorted({canon((0,) + r) for r in itertools.combinations_with_replacement(range(n), W - 1)})
# Burnside check of the class count: group k -> 2^s k + c acting on W-multisets of Z_n
def burnside():
    tot = 0; G = 0
    for s in range(m):
        for c in range(n):
            G += 1
            perm = [((k << s) + c) % n for k in range(n)]
            # count multisets of size W fixed by perm: product over cycles, generating function
            seen = [False] * n; cyc = []
            for k in range(n):
                if not seen[k]:
                    L = 0; x = k
                    while not seen[x]: seen[x] = True; x = perm[x]; L += 1
                    cyc.append(L)
            poly = [1] + [0] * W
            for L in cyc:
                new = [0] * (W + 1)
                for d in range(W + 1):
                    if poly[d]:
                        for t in range(0, W + 1 - d, L): new[d + t] += poly[d]
                poly = new
            tot += poly[W]
    assert tot % G == 0
    return tot // G
print(f"classes of {W}-multisets: {len(classes)} (Burnside {burnside()})", flush=True)

ninst = 0; npass = 0; results = []; post = set()
for ks in classes:
    r = sum(PAR[k] for k in ks)
    frs = [FR[k] for k in ks]
    done = set()
    for ch in itertools.product(*[range(len(f)) for f in frs]):
        c = np.zeros(n, dtype=np.int16)
        for l in range(W):
            for z in frs[l][ch[l]]: c[z - 1] += 1
        key = c.tobytes()
        if key in done: continue
        done.add(key); ninst += 1
        D = EV @ c
        if (D < r).any():
            print("PLAIN FAILS DELETION?!", ks, ch); continue
        tight = (D == r) & (r > 0)
        bad = EV[tight].any(axis=0)
        for z in range(1, n + 1):
            if c[z - 1] > 0 and not bad[z - 1]:
                npass += 1
                moves = [w for w in range(1, n + 1) for _ in range(int(c[w - 1]) - (w == z))]
                results.append((ks, ch, z, int(c[z - 1]), int(tight.sum())))
                post.add((ks, tuple(moves)))
                print(f"PASS ks={ks} frames={ch} delete dir {z} (multiplicity {int(c[z-1])}) tight-before {int(tight.sum())} moves={moves}", flush=True)
print(f"block {m} W={W}: classes {len(classes)}, distinct plain multisets {ninst}, passing single deletions {npass}, "
      f"distinct resulting (class, multiset) pairs {len(post)}", flush=True)
