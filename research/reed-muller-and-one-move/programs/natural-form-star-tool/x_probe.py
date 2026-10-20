from xbase import *
import collections
m=5; Fd=Field(m); fo=Forms(m)
sets = {"8a":(0,2,4,8,18,24,25,27), "9a":(0,1,2,5,11,26,27,29,30), "10a":(0,3,6,8,9,18,20,22,28,30)}
for name, ks in sets.items():
    sy = System(fo, canon_targets(Fd, ks))
    print(name, "viol", len(sy.violations()), "tight", len(sy.tight()))
    # slack histogram by weight, for forms with r>0
    h = collections.Counter()
    for w in (8,12,16):
        sel = (fo.wt==w)
        h = collections.Counter(zip(sy.r[sel].tolist()))
        print("  weight",w, "r-hist", sorted(h.items()))
    even = np.nonzero((sy.r==0)&(fo.wt>0))[0]
    print("  forms with r=0:", len(even), "weights", sorted(collections.Counter(fo.wt[even].tolist()).items()))
sy = System(fo, [identity(5)]*8); print("slice 8 in 31: viol", len(sy.violations()), "tight", len(sy.tight()))
sy = System(fo, [identity(5)]*7); print("slice 7 in 31: viol", len(sy.violations()), "tight", len(sy.tight()))
