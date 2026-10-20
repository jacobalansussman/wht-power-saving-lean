"""famI22.py -- w3-asym: famI on cores in the special position Lemma 5 forces: lambda(K p_l^-1 p_j) = 2 for every
ordered pair (rejection sampling; about 6 % of random merges at W = 2).  Same closure-guided frame-matched injection
search on the FULL closure residual.  Usage: famI22.py W nskel seed maxinj nsamp outprefix"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famI
from famI import *

def core22(rng, W, tries=3000):
    for _ in range(tries):
        bs = [prog_letter(rng, [G.ID]) for l in range(W)]
        hs = [G.intern(G.mul(G.D1, b)) for b in bs]
        as_ = [prog_letter(rng, [h]) for h in hs]
        ps = [G.mul(G.D1, a, G.MATS[h]) for a, h in zip(as_, hs)]
        if all(lam(G.mul(G.K, np.conj(ps[l]).T, ps[j])) == 2 for l in range(W) for j in range(W) if l != j):
            break
    else:
        return None
    pids = [G.intern(p) for p in ps]
    rec = [[(l, bs[l])] for l in range(W)]
    rec.append([(W + l, as_[l]) for l in range(W)])
    u = 2 * W
    for l in range(W):
        path = midpaths(rng, pids[l], nsamp=400)
        if not path: return None
        ell, m = path[0]
        rec.append([(u, ell)])
        rec.append([(u + 1 + 2 * l, m)])
    return rec

if __name__ == '__main__':
    famI.run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6], core22)
