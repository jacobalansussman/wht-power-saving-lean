# c9_subcal.py -- calibration of the SUBGAME symmetry screen (letters from B'):
# the plain programs live in the subgame, so every group must let plain 5 in 15 / 6 in 18 / 4 in 12 through in slack mode
# (some array structure); 5 in 13 and 2 in 5 must fail. Then 5 in 14 tight for the small groups (orders 2, 3, 4).
import numpy as np, time
from lib514 import *
import c7_symscreen2 as S7
import c8_s4 as C8
t0 = time.time()
Lsubs = S7.subgroups_upto_conj(LEVI)
for H in Lsubs:
    if len(H) in (1, 7, 21, 168): continue
    info = C8.analyse_subgame(H, S7.Y)
    pts, planes = C8.fixed_subspaces(H)
    cal = {}
    for (W, M) in [(4, 12), (5, 15), (6, 18), (8, 24), (12, 36), (24, 72), (5, 13), (2, 5)]:
        res, lines = S7.describe(info, W, M, tight=False, gauge=False)
        cal[f'{W}in{M}'] = len(res)
    res, lines = S7.describe(info, 5, 14, tight=True, gauge=False)
    print(f'H order {len(H)} fixes z-vectors {pts} z-planes {planes}: slack-mode array structures passing {cal}; tight 5in14 shapes {sum(len(s) for _,_,s in res)}', flush=True)
    for l in lines[:4]: print('      ', l)
print('[%.0fs]' % (time.time() - t0))
