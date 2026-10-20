"""vsp.py -- v-o3-close: an INDEPENDENT symplectic model of the block-4 engine (no code from sp.py).

Vectors of F_2^8 are ints v = (x << 4) | z  (x = X-part, address bits 0..3; z = Z-part).  Form <v,w> = x.z' + z.x'.
Lagrangians are enumerated by BRUTE FORCE (all 4-dim totally isotropic subspaces), not as an orbit.
Model (to be checked against matrices in vmat.py):
  sigma(D_{e1}) = transvection by a = (e1|0);  sigma(K) = product of the transvections by (e_j|0);
  depth(g) = dist(sigma(g) L0, L0);  dist_K(g) = dist(sigma(g) Lambda_K, L0),  Lambda_K = sigma(K)^-1 L0;
  dist(A, B) = 4 - dim(A cap B)."""
import itertools

def X(v): return v >> 4
def Z(v): return v & 15
def pc(a): return bin(a).count('1')
def form(v, w):
    return (pc(X(v) & Z(w)) + pc(Z(v) & X(w))) & 1

def spanset(vs):
    S = {0}
    for v in vs:
        if v not in S:
            S = S | {s ^ v for s in S}
    return frozenset(S)

def dim(S):
    return len(S).bit_length() - 1

def dist(A, B):
    return 4 - dim(A & B)

ALLV = list(range(1, 256))

def all_lagrangians():
    """brute force: grow totally isotropic subspaces one vector at a time"""
    layer = {frozenset([0])}
    for k in range(4):
        nxt = set()
        for S in layer:
            perp = [v for v in ALLV if v not in S and all(form(v, s) == 0 for s in S)]
            for v in perp:
                nxt.add(frozenset(S | {s ^ v for s in S}))
        layer = nxt
    return layer

def transvection(a):
    return lambda v: v ^ a if form(v, a) else v

def apply(f, S):
    return frozenset(f(v) for v in S)

L0 = spanset([1 << j for j in range(4)])                 # Z-type: x = 0
A1 = 1 << 4                                              # (e1 | 0)
tau = transvection(A1)
def sigK(v):
    for j in range(4):
        v = transvection(1 << (4 + j))(v)
    return v
LK = apply(sigK, L0)            # sigma(K) is an involution, so sigma(K)^-1 L0 = sigma(K) L0
LD = apply(tau, L0)

if __name__ == '__main__':
    LAG = all_lagrangians()
    from collections import Counter
    print("Lagrangians (brute force):", len(LAG))
    print("distance profile from L0:", sorted(Counter(dist(L, L0) for L in LAG).items()))
    print("dist(L0, Lambda_K) =", dist(L0, LK), "; dist(L0, Lambda_D) =", dist(L0, LD))
    print("Lambda_K == {(z|z)}:", LK == frozenset((z << 4) | z for z in range(16)))
