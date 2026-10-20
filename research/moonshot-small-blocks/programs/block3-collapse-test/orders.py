"""orders.py -- ms-xw3.  Orders (as permutations, ignoring phases) of a few XW free words at f = 2..6: a word g with
order growing exponentially in f makes 1 + eps g a bounded mix whose inverse sum_k (-eps)^k g^k is NOT bounded."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from math import gcd
from xw3lib import Game


def perm_order(sig):
    n = sig.size; seen = np.zeros(n, bool); L = 1
    for s in range(n):
        if seen[s]: continue
        c = 0; x = s
        while not seen[x]:
            seen[x] = True; x = sig[x]; c += 1
        L = L * c // gcd(L, c)
    return L


words = ["A", "a0", "E01", "T2", "A.T2", "T2.A", "a0.T2", "E01.T2", "A.T2.B", "A.C01", "A.P01", "E01.E12", "A.T2.P01", "a0.a1.T2"]
G = Game([2, 3, 4, 5, 6], W=1, S=1, seed=1)
print("word".ljust(12), " ".join(f"f={f}".rjust(10) for f in G.FS))
for w in words:
    g = G.word(w)
    print(w.ljust(12), " ".join(str(perm_order(g[f].sig)).rjust(10) for f in G.FS), flush=True)
