import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from asym import *
P = plain(2)
print(P.describe())
DT = DTest(P)
print("plain", DT.residual(np.zeros(0), detail=True))
print("plain full2", [r[0] for r in full2(P, np.zeros(0))])
