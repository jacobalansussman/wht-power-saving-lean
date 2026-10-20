# Which sets S of nonzero directions z in F_2^m have sum_z z z^T = 0 (mod 2)?  e7-rm's L-det/L-out need the
# smallest such "dependency supports"; claim: in F_2^4 the smallest have 7 elements and are exactly the 15 3-flats
# (nonzero points of 3-dim subspaces).  Exhaustive over all subsets of size <= 7 (block 4) and <= 7 (block 5).
import itertools, sys
def zzT(z, m):
    v = 0; k = 0
    for i in range(m):
        for j in range(i, m):
            if (z >> i) & 1 and (z >> j) & 1: v |= 1 << k
            k += 1
    return v
def is_flat(S, m):
    s = set(S) | {0}
    return all((a ^ b) in s for a in s for b in s) and len(s) & (len(s) - 1) == 0
for m in (3, 4, 5):
    Z = list(range(1, 1 << m)); T = {z: zzT(z, m) for z in Z}
    for size in range(1, 8):
        hits = []
        for S in itertools.combinations(Z, size):
            acc = 0
            for z in S: acc ^= T[z]
            if acc == 0: hits.append(S)
        if hits:
            nflat = sum(is_flat(S, m) for S in hits)
            print(f"block {m}: size {size}: {len(hits)} dependency supports, {nflat} of them are flats (nonzero points of a subspace)")
    print(f"block {m}: done", flush=True)
