# hand_b4.py: search for a HAND-READABLE certificate of Theorem 17. For each of the 7 affine planes P of the even-weight
# subspace E (P avoids 0), average each piece's slack bound uniformly over ALL odd planes containing the piece, give the
# planes of the two S_4-orbits weights lamA, lamB, and compute the coefficient of every c_z. If all coefficients are equal
# (= kappa), summing gives  kappa * sum_z c_z >= 16 (7 terms weighted)  =>  sum_z c_z >= 4.
import itertools
from fractions import Fraction as Fr
pts = range(1, 16)
def is_plane(S): S = sorted(S); return len(S) == 4 and S[0] ^ S[1] ^ S[2] ^ S[3] == 0
planes = [frozenset(S) for S in itertools.combinations(pts, 4) if is_plane(S)]
basis = {1, 2, 4, 8}
odd_planes = [P for P in planes if len(P & basis) % 2 == 1]
Estar = [z for z in pts if bin(z).count("1") % 2 == 0]
Eplanes = [P for P in planes if P <= set(Estar)]
print("odd planes", len(odd_planes), "even planes of E", len(Eplanes), [sorted(P) for P in Eplanes])
def contrib(P):
    co = {z: Fr(0) for z in pts}
    for z in P: co[z] -= 1                      # -c(P)
    for x in P:                                  # diagonal pieces: average over odd planes containing x
        G = [Q for Q in odd_planes if x in Q]
        for Q in G:
            for z in Q: co[z] += Fr(1, len(G))
    for a, b in itertools.combinations(sorted(P), 2):
        G = [Q for Q in odd_planes if a in Q and b in Q]
        for Q in G:
            for z in Q: co[z] += Fr(2, len(G))
    return co
orbA = [P for P in Eplanes if 15 in P]; orbB = [P for P in Eplanes if 15 not in P]
cA = {z: sum(contrib(P)[z] for P in orbA) for z in pts}; cB = {z: sum(contrib(P)[z] for P in orbB) for z in pts}
print("orbit A (contain 15):", len(orbA), " orbit B:", len(orbB))
print("coefficients, orbit A sum:", {z: str(v) for z, v in cA.items()})
print("coefficients, orbit B sum:", {z: str(v) for z, v in cB.items()})
# the hand certificate: weight 18 on the 4 planes containing 15 = 1111, weight 17 on the other 3
lamA, lamB = 18, 17
tot = {z: lamA * cA[z] + lamB * cB[z] for z in pts}
const = 16 * (lamA * len(orbA) + lamB * len(orbB))      # each plane-term: 4 diagonal pieces (-1 each) + 6 pairs x 2 (-1 each)
print("combined coefficients:", sorted(set(str(v) for v in tot.values())), " constant:", const)
k = list(tot.values())[0]
if all(v == k for v in tot.values()):
    print(f"HAND CERTIFICATE: {k} * sum_z c_z >= {const}, i.e. sum_z c_z >= {Fr(const, 1) / k}")
