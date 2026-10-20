import numpy as np
from vtests import Tester
t = Tester(5)
sets = {"rm-theory 8": (0,2,4,8,18,24,25,27), "additive 8": (0,0,1,9,16,18,21,25), "Frobenius 8": (0,0,0,1,2,4,8,16),
        "N1": (0,8,9,10,14,17,23,28), "rm-theory 9": (0,1,2,5,11,26,27,29,30), "C1uC3 10": (1,2,4,8,16,3,6,12,17,24),
        "slice 8": (0,)*8, "slice 7": (0,)*7}
E = t.E.astype(np.int64); D = E.sum(1)
for name, ks in sets.items():
    r = sum(t.parity(k).astype(np.int64) for k in ks)
    print(f"{name:12s} one move per direction: violations {int((D < r).sum())}, tight forms {int(((D == r) & (r > 0)).sum())}")
