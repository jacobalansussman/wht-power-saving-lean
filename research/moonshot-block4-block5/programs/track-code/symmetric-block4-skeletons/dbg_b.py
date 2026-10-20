import sys
import symsolve as S
from symsolve import *
orig = Skel.add_move
def am(self, kind, summ, Dm_override=None):
    r = orig(self, kind, summ, Dm_override)
    print("  added", kind, "summands", len(summ), "orbit", r, "terms", self.moves[-1]["n"])
    return r
Skel.add_move = am
for seed in range(3):
    sk, cost = build_skeleton(3, 2, 'Z2swap', sys.argv[1].split(','), 3, 700 + seed)
    print("seed", seed, "->", cost)
