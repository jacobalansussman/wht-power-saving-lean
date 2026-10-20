#!/usr/bin/env python3
"""fc-reopt: check the suffix-interval capacity condition on an integral class solution (SOL-*.json of build_fc.py):
for I = [A, N-1] and each set S of widths, sum over bins of min(slots of S, |I n window|) >= (roles of S per call) |I|.
usage: diag_end.py LINE(3945|2979) HUB3815 ORDER SOL.json [LAST=80]"""
import sys, os, json, gzip
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hall_fc as HF
from collections import Counter
line, hub, orderp, solp = sys.argv[1:5]; LAST = int(sys.argv[5]) if len(sys.argv) > 5 else 80
order = json.loads(open(orderp).read()); T = 1 + max(r for _, r in order)
if line == '3945':
    sup = [26, 2, 18, 53, 96, 213]; local = {1: 2440, 2: 725, 3: 224, 4: 102, 5: 20, 8: 8}
else:
    prov = json.load(open(os.path.join(os.path.dirname(orderp), 'PROVIDERS.json')))
    inv = Counter(prov['gauge44']['widths']); sup = [inv.get(w, 0) for w in (4, 5, 6, 7, 8, 9)]
    local = {int(k): v for k, v in prov['gauge44']['actual_paid_profile'].items()}
H = HF.load(os.path.join(hub, 'research/staged544-clock-ceiling/construction/hall_lp.py'), T, sup, local, 24, False, line == '3945')
cap, lo, hi = H.classes(order)
N = H.N
sol = json.load(open(solp))
cols = [(int(k), int(j), float(v)) for k, j, v in sol['cols'] if v > 1e-9]
CK = np.array([k for k, j, v in cols]); CN = np.array([v for k, j, v in cols])
CS = np.array([[sum(int(H.PATS[int(cap[k])][0][j][w]) for w in range(6) if S >> w & 1) for S in H.SUBSETS] for k, j, v in cols])
clo = lo[CK] + 1; chi = np.minimum(hi[CK], N - 1)
viol = []
for A in range(N - LAST, N):
    ov = np.maximum(chi - np.maximum(clo, A) + 1, 0)
    capS = (CN[:, None] * np.minimum(CS, ov[:, None])).sum(0)
    need = H.SUBSET_SUP * (N - A)
    badi = np.where(capS < need - 1e-9)[0]
    if len(badi):
        viol.append((A, [(int(H.SUBSETS[i]), float(capS[i]), int(need[i])) for i in badi[:4]]))
print('suffix starts violating the interval condition (A, [(subset mask over widths 4..9, capacity, need)]):')
for v in viol[:20]: print('  ', v)
print('total violating starts:', len(viol), 'of', LAST)
