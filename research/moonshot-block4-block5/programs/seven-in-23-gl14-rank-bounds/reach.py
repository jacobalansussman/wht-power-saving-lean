"""reach.py -- w8-bpgl14: which of the 30 slots (neighbours of F) can terms OTHER than the top's F1/F2-part reach in output 0?
Letters keep cells.  Output 0 reads: the inputs (cell L_Z), O1 (lam1), O2 (t1), F1, F2 (F) and the top, whose own other parts
are D3 b w for w = inputs, O1, O2: cells = neighbours of L_Z, lam1, t1.  Prints the slot sets."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4test import *
L = Lift(seed=21, psi=1); D = D4(L); P = L.P
sl = D.slot
dist = lambda d: int(P.dZ[sl[d]])
nb = lambda c: [d for d in range(30) if P.D[sl[d], c] == 1]
R_t1 = set(nb(L.t1)) | set(d for d in range(30) if sl[d] == L.t1)
R_l1 = set(nb(L.lam1)); R_Z = set(nb(L.iZ))
print('slot distances to the input cell:', [dist(d) for d in range(30)])
print('S (distance 4):', sorted(D.S), ' d0', D.d0)
print('d4-free (not distance 4):', sorted(d for d in range(30) if d not in D.S))
print('reached via t1 (O2 itself or the top\'s O2-part):', sorted(R_t1))
print('reached via lam1 (the top\'s O1-part):', sorted(R_l1), ' via L_Z:', sorted(R_Z))
print('union:', sorted(R_t1 | R_l1 | R_Z))
