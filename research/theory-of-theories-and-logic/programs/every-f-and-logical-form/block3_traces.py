"""block3_traces.py -- v-ml-logic: sample traces of the block-3 engine group (own generator code via
engine_group_exact.gens, converted to numpy) with random words of length 1..80; exact traces lie in Z[i,1/2] and are
rounded to 1/1024. Prints the distinct nonzero traces and |tr|^2 values (ml-logic claims r'(3) >= 28; the Clifford
trace fact would give exactly 4(2m+1) = 28).  Usage: python3 -I block3_traces.py  (tiny, < 1 min)"""
import sys, random
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import engine_group_exact as E
m = 3
gs = [np.array([[complex(*e) for e in row] for row in N]) / 2**k for (N, k) in E.gens(m).values()]
rng = random.Random(11); vals = set()
for t in range(40000):
    g = np.eye(8, dtype=complex)
    for _ in range(rng.randrange(1, 81)): g = gs[rng.randrange(len(gs))] @ g
    tr = np.trace(g); vals.add((round(tr.real * 1024), round(tr.imag * 1024)))
nz = [v for v in vals if v != (0, 0)]
mods = sorted({(a * a + b * b) / 1024**2 for a, b in vals})
print(f"m = 3: {len(vals)} distinct traces sampled, {len(nz)} nonzero; |tr|^2 values {mods}")
