"""frames of every copy target (5-sets of directions whose forms add up to the target's form), my code."""
import itertools
from vcore import EXP, dir_vec, tgt_vec
DV = {x: dir_vec(x) for x in range(1, 32)}
TV = {k: tgt_vec(EXP[(2 * k) % 31]) for k in range(31)}
inv = {v: k for k, v in TV.items()}
frames = {k: [] for k in range(31)}
for S in itertools.combinations(range(1, 32), 5):
    v = 0
    for x in S: v ^= DV[x]
    if v in inv: frames[inv[v]].append(S)
