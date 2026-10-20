"""selective_sp.py -- o3-close: EXHAUSTIVE test of o3-build's "no selective advance" observation in the exact
symplectic model (sp.py, cross-checked against the matrices by check_model.py).

Joint move u3 = D(r u1 + r' u2), u1 = D b1 y1, u2 = D b2 y2 (single letters).  Its depth-2 terms h1 (y1), h2 (y2) have
states (mu_j, kappa_j) = tau sigma(r) tau sigma(b1) (L0, Lambda_K) etc.  A later letter t acts on both through
N = sigma(t)^-1 Lambda_D (a neighbour of L0, 30 choices):  depth(D t h) = dist(mu, N), dist_K(D t h) = dist(kappa, N).
State of h after r, b1: (tau N1, tau Y1) with N1 in nbr(L0), Y1 opposite N1 (all such pairs occur).
SELECTIVE (o3-build's definition): D t h1 at depth 3 next to K (dist(mu1,N)=3, dist(kappa1,N)=1) while D t h2 has
depth <= 2; and the depth-only version (D t h1 depth 3, D t h2 depth <= 2)."""
import sys, os, itertools
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sp
from sp import dist, L0, LD, img, TAU

LAG = sorted(sp.all_lagrangians(), key=lambda S: sorted(S))
NB = [L for L in LAG if dist(L, L0) == 1]
OPP = {}
def opp(X):
    k = X
    if k not in OPP:
        OPP[k] = [Y for Y in LAG if dist(X, Y) == 4]
    return OPP[k]
tau = lambda S: img(TAU, S)

# depth-2 single-letter terms: mu = tau N1 with dist(N1, LD) = 2
N2 = [N for N in NB if dist(N, LD) == 2]
print("neighbours N1 of L0 with tau N1 at depth 2:", len(N2), "(of 30)")
states = []                    # (mu, kappa) of depth-2 terms D r D b
for N1 in N2:
    for Y in opp(N1):
        states.append((tau(N1), tau(Y)))
print("depth-2 single-letter term states:", len(states))
geo_states = [s for s in states if dist(s[1], L0) == 2]
print("  of which on a geodesic to K (dist_K = 2):", len(geo_states))
mus2 = sorted({tau(N1) for N1 in N2}, key=lambda S: sorted(S))
print("distinct mu for depth-2 terms:", len(mus2))

def I_of(mu):
    return mu & L0

stat = Counter(); stat_depth = Counter(); crit = Counter()
for (mu1, ka1) in geo_states:
    for mu2 in mus2:
        sel = [N for N in NB if dist(mu1, N) == 3 and dist(ka1, N) == 1 and dist(mu2, N) <= 2]
        seld = [N for N in NB if dist(mu1, N) == 3 and dist(mu2, N) <= 2]
        stat[len(sel)] += 1
        stat_depth[len(seld)] += 1
        J1 = ka1 & L0; I2 = I_of(mu2); I1 = I_of(mu1)
        crit[(len(sel) > 0, sp.dim(J1 & I2), I1 == I2)] += 1
print("pairs (geodesic term h1, depth-2 term h2):", sum(stat.values()))
print("  number of letter classes N giving a SELECTIVE GEODESIC advance of h1:", sorted(stat.items()))
print("  number giving a selective DEPTH advance of h1 (any frame at depth 3):", sorted(stat_depth.items()))
print("  (selective exists, dim(J1 cap I2), I1 == I2) :", sorted(crit.items()))
# check the claimed characterisations on all 30 N and all states
ok1 = ok2 = 0
for (mu1, ka1) in states:
    I1 = mu1 & L0
    for N in NB:
        H = N & L0
        if (dist(mu1, N) <= 2) == (I1 <= H):
            ok1 += 1
for (mu1, ka1) in geo_states:
    J1 = ka1 & L0
    closer = [N for N in NB if dist(ka1, N) == 1]
    if sorted(sorted(N & L0) for N in closer) == sorted(sorted(H) for H in {N & L0 for N in NB} if J1 <= H):
        ok2 += 1
print(f"check 'h not advanced <=> N cap L0 contains mu cap L0': {ok1} of {len(states) * 30}")
print(f"check 'the 3 letters N advancing a geodesic term towards K are those with N cap L0 = the 3 hyperplanes "
      f"containing J = kappa cap L0' (as sets of hyperplanes): {ok2} of {len(geo_states)}")
