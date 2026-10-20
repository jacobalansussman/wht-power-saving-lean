"""vscan_alt.py: count plain multisets and passing sets per class representative WITHOUT the stabiliser reduction
(distinct sorted multisets per target-class representative), to compare with the track's 7,076."""
import itertools
import numpy as np
exec(open("vscan.py").read().split("ms = sorted")[0])
reps = sorted(set(canon(c, [])[0] for c in itertools.combinations_with_replacement(range(31), 3)))
np_ = 0; npass = 0; full = set(); fullp = set()
for ks in reps:
    seen = set(); seenp = set()
    for fr in itertools.product(*(range(len(frames[k])) for k in ks)):
        moves = sorted(x for k, f in zip(ks, fr) for x in frames[k][f])
        seen.add(tuple(moves))
        cnt = sum(DT[x].astype(np.int32) for x in moves); r = sum(TT[k].astype(np.int32) for k in ks)
        tight = (cnt - r) == 0
        for d in sorted(set(moves)):
            if (tight & (DT[d] == 1)).any(): continue
            m14 = list(moves); m14.remove(d); seenp.add(tuple(m14))
    np_ += len(seen); npass += len(seenp)
print("per-rep distinct plain multisets:", np_, " per-rep distinct passing 14-multisets:", npass)
