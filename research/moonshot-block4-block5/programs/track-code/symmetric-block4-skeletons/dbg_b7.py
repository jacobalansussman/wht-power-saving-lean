import sys
from symsolve import *
orig = Skel.add_move
def am(self, kind, summ, Dm_override=None):
    r = orig(self, kind, summ, Dm_override)
    d = Dm_override
    fixd = None
    if d is not None:
        fixd = all(self.G.is_monomial(self.G.mul(c, d, np.conj(c).T, np.conj(d).T)) for (c, pi, cp) in self.R.Hc[1:])
    print("  added", kind, "summands", len(summ), "orbit", r, "terms", self.moves[-1]["n"], "D fixed:", fixd)
    return r
Skel.add_move = am
sk, cost = build_skeleton(4, int(sys.argv[2]), sys.argv[3], sys.argv[1].split(','), 1, 4100)
print("cost", cost)
