"""vcensus.py -- v-o3-close: independent re-derivation of o3-close Proposition S's census.

P0 generators are sigma(b) of random 16x16 monomial Cliffords (vmat.py), so the state set below is built from the
real letters, not from the description "kappa = any opposite of mu".
States of h = D r D b (single letters): (mu, kappa) = sigma(h)(L0, Lambda_K). Built as
  S_Db   = { tau p (L0, Lambda_K) : p in P0 }              (P0-orbit of (L0, Lambda_K), then tau)
  S_DrDb = { tau p s : p in P0, s in S_Db }                 (P0-orbit of S_Db, then tau)
A letter t (class N = sigma(t)^-1 Lambda_D, a neighbour of L0) sends state (mu, kappa) to depth dist(mu, N) and
dist_K dist(kappa, N) (checked at matrix level in vselect.py).
Census: geodesic h1 (depth 2, dist_K 2) x depth-2 forward frames mu2; count classes N with dist(mu1,N)=3,
dist(kappa1,N)=1, dist(mu2,N)<=2.  Also the joint distribution with dim(J1 cap I2), and checks (a), (b)."""
import sys, os, random, time
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vsp, vmat
from vsp import dist, dim, L0, LK, LD
from collections import Counter

t0 = time.time()
rng = random.Random(5)
GENS = [vmat.sigma(vmat.rand_free(rng)) for _ in range(6)]
def act(tab, S): return frozenset(tab[v] for v in S)
TAU = [vsp.tau(v) for v in range(256)]

def orbit_pairs(seeds):
    seen = set(seeds); fr = list(seeds)
    while fr:
        nf = []
        for (A, B) in fr:
            for g in GENS:
                p = (act(g, A), act(g, B))
                if p not in seen:
                    seen.add(p); nf.append(p)
        fr = nf
    return seen

O1 = orbit_pairs([(L0, LK)])
print("P0-orbit of (L0, Lambda_K):", len(O1), "(expect 1024 opposites of L0)")
S_Db = {(act(TAU, A), act(TAU, B)) for (A, B) in O1}
O2 = orbit_pairs(list(S_Db))
S_h = {(act(TAU, A), act(TAU, B)) for (A, B) in O2}
print("states of D r D b:", len(S_h), " by depth:", sorted(Counter(dist(m, L0) for (m, k) in S_h).items()))
dep2 = [(m, k) for (m, k) in S_h if dist(m, L0) == 2]
geo = [(m, k) for (m, k) in dep2 if dist(k, L0) == 2]
mus2 = sorted({m for (m, k) in dep2}, key=lambda S: sorted(S))
NB = sorted({act(g, LD) for g in GENS} | set(), key=lambda S: sorted(S))
# all 30 neighbours: P0-orbit of Lambda_D
seen = {LD}; fr = [LD]
while fr:
    nf = []
    for A in fr:
        for g in GENS:
            B = act(g, A)
            if B not in seen:
                seen.add(B); nf.append(B)
    fr = nf
NB = sorted(seen, key=lambda S: sorted(S))
print(f"depth-2 states {len(dep2)}, geodesic {len(geo)}, distinct mu2 {len(mus2)}, letter classes (neighbours) {len(NB)}")
assert all(dist(N, L0) == 1 for N in NB)

nsel = Counter(); joint = Counter(); both = Counter()
for (mu1, ka1) in geo:
    J1 = ka1 & L0; I1 = mu1 & L0
    for mu2 in mus2:
        I2 = mu2 & L0
        sel = [N for N in NB if dist(mu1, N) == 3 and dist(ka1, N) == 1 and dist(mu2, N) <= 2]
        nsel[len(sel)] += 1
        joint[(len(sel), dim(J1 & I2))] += 1
print("pairs:", sum(nsel.values()))
print("number of selective geodesic classes per pair:", sorted(nsel.items()))
print("(number of classes, dim(J1 cap I2)):", sorted(joint.items()))
# (a): dist(mu, N) <= 2  <=>  N cap L0 contains mu cap L0   (all depth-2 states, all N)
ok = sum((dist(m, N) <= 2) == ((m & L0) <= (N & L0)) for (m, k) in dep2 for N in NB)
print(f"(a) holds in {ok} of {len(dep2) * len(NB)}")
# (b): N nearer K (dist(kappa, N) = 1) <=> N cap L0 contains J and N is the unique nearer point of its line
okb = 0
for (m, k) in geo:
    J = k & L0
    nearer = [N for N in NB if dist(k, N) == 1]
    Hs = sorted(sorted(N & L0) for N in nearer)
    cont = sorted(sorted(H) for H in {N & L0 for N in NB} if J <= H)
    okb += (len(nearer) == 3 and Hs == cont)
print(f"(b) holds in {okb} of {len(geo)}")
# the BOTH-WAYS condition (needed by shape R): pairs of geodesic h1, h2 that come from ONE u3 (independent r, r')
# fraction with selective letters in both directions, over all (state1, state2) pairs with both geodesic
cnt = Counter()
geo_by_mu = {}
for (m, k) in geo:
    geo_by_mu.setdefault(m, []).append(k)
smp = random.Random(9)
for _ in range(200000):
    m1, k1 = smp.choice(geo); m2, k2 = smp.choice(geo)
    d1 = dim((k1 & L0) & (m2 & L0)) > 0; d2 = dim((k2 & L0) & (m1 & L0)) > 0
    cnt[(d1, d2)] += 1
print("random geodesic pairs (h1, h2), selective (1->, 2->):", sorted(cnt.items()))
print(f"time {time.time() - t0:.1f}s")
