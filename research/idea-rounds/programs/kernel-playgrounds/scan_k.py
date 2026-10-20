import sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hookk import hook
P1, P2 = 2147483647, 2147483629
K = {'H': [[1,1],[1,-1]], 'L': [[1,0],[1,1]], 'D7': [[2,3],[5,7]]}
m = int(sys.argv[1])
for k in [int(a) for a in sys.argv[2:]]:
    for name, kk in K.items():
        r = hook(kk, m, k, P1); r2 = hook(kk, m, k, P2)
        print(f"m={m} k={k} f={k+1} dimN={r['dimN']} {name}: cap={r['cap']} dem={r['dem']} ratio={r['dem']/r['cap']:.4f} agree={r==r2}", flush=True)
