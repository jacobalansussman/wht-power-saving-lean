"""packing.py: for a canonical target set at block 5, the largest number of arrays whose frames are pairwise
disjoint (the plain packing number in the one-move-per-direction game), and the disjointness graph."""
import sys, itertools
from inst import Field
fd = Field(5)
def frames(ks): return [set(fd.frame(k)) for k in ks]
def packing(ks):
    fr = frames(ks); W = len(ks); best = []
    for r in range(W, 0, -1):
        for S in itertools.combinations(range(W), r):
            if all(not (fr[i] & fr[j]) for i, j in itertools.combinations(S, 2)):
                return r, S
    return 0, ()
if __name__ == "__main__":
    for spec in sys.argv[1:]:
        ks = [int(v) for v in spec.split(",")]
        print(spec, "packing number, a maximum family:", packing(ks), flush=True)
