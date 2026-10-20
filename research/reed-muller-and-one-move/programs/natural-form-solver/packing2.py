"""packing2.py: ALL frames of each canonical target (5-sets of directions z with sum z z^T = T_k over F2) and the
true plain packing number of a target multiset in the one-move-per-direction game (pairwise disjoint frames)."""
import sys, os, itertools
from inst import Field, mons
MB = int(os.environ.get("HM", "5"))
fd = Field(MB)
ML = mons(MB)
def smask_dir(z):
    b = [(z >> j) & 1 for j in range(MB)]
    return sum((b[i] & b[j]) << q for q, (i, j) in enumerate(ML))
def smask_T(T):
    return sum((T[i][j] & 1) << q for q, (i, j) in enumerate(ML))
DM = {z: smask_dir(z) for z in range(1, 1 << MB)}
bucket = {}
for S in itertools.combinations(range(1, 1 << MB), MB):
    v = 0
    for z in S: v ^= DM[z]
    bucket.setdefault(v, []).append(frozenset(S))
def all_frames(k):
    return bucket.get(smask_T(fd.target(k)), [])
def packing(ks):
    fr = [all_frames(k) for k in ks]; W = len(ks); best = [0, None]
    def rec(i, used, cnt, choice):
        if cnt + (W - i) <= best[0]: return
        if i == W:
            best[0] = cnt; best[1] = list(choice); return
        for F in fr[i]:
            if not (F & used): rec(i + 1, used | F, cnt + 1, choice + [(i, sorted(F))])
        rec(i + 1, used, cnt, choice)
    rec(0, frozenset(), 0, [])
    return best
if __name__ == "__main__":
    print("frames per target:", {k: len(all_frames(k)) for k in range(3)})
    for spec in sys.argv[1:]:
        ks = [int(v) for v in spec.split(",")]
        b = packing(ks)
        print(spec, "true plain packing number", b[0], "e.g.", b[1], flush=True)
