"""famS.py -- w3-asym: the case Theorem R does NOT cover: a sigma IN THE MERGE.  Core: starts s_l = D b_l x_l; merge
u = D(a_0 s_0 + c a_0 sigma s_0 + sum_(l>=1) a_l s_l) with sigma an order-3 address map (a cyclic shift of three address
bits, conjugated by a random GL(4,2) element) chosen so that both array-0 merge terms progress (layer 2); then the
private paths of famI (for array 0 along the first merge term) and the same closure-guided frame-matched injection
search on the FULL closure residual.  Usage: famS.py W nskel seed maxinj nsamp outprefix"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famI
from famI import *

def order3_letters(rng, n3=40):
    """order-3 address maps: conjugates of the 3-cycle of bits (0 1 2) by random GL(4,2) address maps"""
    base = [0] * 16
    for x in range(16):
        b = [(x >> j) & 1 for j in range(4)]
        y = b[2] | (b[0] << 1) | (b[1] << 2) | (b[3] << 3)
        base[x] = y
    S = G.mono(base, [0] * 16)
    out = [S]
    for _ in range(n3):
        while True:
            A = [rng.randrange(1, 16) for _ in range(4)]
            span = {0}; ok = True
            for a in A:
                if a in span: ok = False; break
                span |= {x ^ a for x in span}
            if ok: break
        sig = []
        for x in range(16):
            y = 0
            for j in range(4):
                if (x >> j) & 1: y ^= A[j]
            sig.append(y)
        M = G.mono(sig, [0] * 16)
        out.append(G.mul(M, S, np.conj(M).T))
    return out

def core_sigma(rng, W, tries=300):
    sigs = order3_letters(rng)
    for _ in range(tries):
        bs = [prog_letter(rng, [G.ID]) for l in range(W)]
        hs = [G.intern(G.mul(G.D1, b)) for b in bs]
        as_ = [prog_letter(rng, [h]) for h in hs]
        ps = [G.mul(G.D1, a, G.MATS[h]) for a, h in zip(as_, hs)]
        if not all(lam(G.mul(G.K, np.conj(ps[l]).T, ps[j])) <= 3 for l in range(W) for j in range(W) if l != j):
            continue
        rng.shuffle(sigs)
        sg = None
        for s in sigs:
            q = G.intern(G.mul(G.D1, as_[0], s, G.MATS[hs[0]]))
            if G.dist(q) == 2 and lam(G.MATS[q]) == 2 and G.frm(q) != G.frm(G.intern(ps[0])):
                sg = s; break
        if sg is not None:
            break
    else:
        return None
    pids = [G.intern(p) for p in ps]
    rec = [[(l, bs[l])] for l in range(W)]
    rec.append([(W, as_[0]), (W, G.mul(as_[0], sg))] + [(W + l, as_[l]) for l in range(1, W)])
    u = 2 * W
    for l in range(W):
        path = midpaths(rng, pids[l], nsamp=400)
        if not path: return None
        ell, m = path[0]
        rec.append([(u, ell)])
        rec.append([(u + 1 + 2 * l, m)])
    return rec

if __name__ == '__main__':
    famI.run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6], core_sigma)
