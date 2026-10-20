#!/usr/bin/env python3
"""Independent check of t-sandwich's m=3 witness (no import of frames.py / twist.py).
Frames typed in from out/witness_m3.txt: printed as (x0 x1 x2 | z0 z1 z2).  Convention of frames.py:
L_U = {(u | u+v): u in U, v in U^perp}, L_0={(0|z)}, L_F={(z|z)}.  Ports (odd u): X_u = I(L_<u>, L_F), Y_u = I(L_0, L_{u^perp})."""
import itertools
m = 3; M = 7
def popc(a): return bin(a).count("1")
def om(u, v): return (popc((u & M) & (v >> m)) + popc((v & M) & (u >> m))) & 1
def span(vs):
    s = {0}
    for v in vs:
        if v not in s: s |= {x ^ v for x in s}
    return frozenset(s)
def dim(S): return len(S).bit_length() - 1
def vec(s):   # "(011|000)" printed order x0x1x2|z0z1z2
    xs, zs = s.strip("()").split("|")
    x = sum(int(c) << i for i, c in enumerate(xs)); z = sum(int(c) << i for i, c in enumerate(zs))
    return x | (z << m)
def fr(*bs):
    S = span([vec(b) for b in bs]); assert len(S) == 8 and all(om(a, b) == 0 for a in S for b in S); return S
def perp(U): return frozenset(v for v in range(8) if all(popc(u & v) % 2 == 0 for u in U))
def L(U): U = span(U); return frozenset(u | ((u ^ v) << m) for u in U for v in perp(U))
def lab(s): return sum(int(c) << i for i, c in enumerate(s))   # label in printed order
d = lambda A, B: m - dim(A & B)
L0 = frozenset(z << m for z in range(8)); LF = frozenset(z | (z << m) for z in range(8))
def interval(A, B, allL): return {Z for Z in allL if d(A, Z) + d(Z, B) == d(A, B)}
# all Lagrangians
level = {frozenset([0])}
for k in range(m):
    level = {frozenset(S | {x ^ v for x in S}) for S in level for v in range(1, 64) if v not in S and all(om(v, w) == 0 for w in S)}
allL = level; assert len(allL) == 135
uS, uT = lab("100"), lab("111")      # code's u_S=0b001 is x0, printed "100"
print("u_S.u_T =", popc(uS & uT) % 2)
XS = interval(L([uS]), LF, allL); YS = interval(L0, L(perp([uS])), allL)
XT = interval(L([uT]), LF, allL); YT = interval(L0, L(perp([uT])), allL)
print("|X_S|,|Y_S|,|X_T|,|Y_T| =", len(XS), len(YS), len(XT), len(YT), " X_T & Y_S empty (old: S does not know T):", not (XT & YS))
f57 = fr("(011|000)", "(100|100)", "(000|011)"); f134 = fr("(111|100)", "(111|010)", "(111|001)")
f71 = fr("(000|100)", "(000|010)", "(000|001)"); f47 = fr("(101|000)", "(000|010)", "(000|101)")
f55 = fr("(011|000)", "(000|100)", "(000|011)"); f51 = fr("(101|000)", "(010|010)", "(000|101)")
f79 = fr("(100|100)", "(000|010)", "(000|001)"); f132 = fr("(111|100)", "(101|010)", "(111|001)")
names = {f57: "57", f134: "134", f71: "71", f47: "47", f55: "55", f51: "51", f79: "79", f132: "132"}
ident = {f57: "L_<100,011>", f134: "L_<111>", f71: "L_0", f47: "L_<101>", f55: "L_<011>", f51: "L_<101,010>", f79: "L_<100>"}
for f, s in ident.items():
    U = [lab(t) for t in s[3:-1].split(",")] if s != "L_0" else []
    assert L(U) == f, s
print("frames identified as subspace frames:", ", ".join(f"{names[f]}={s}" for f, s in ident.items()))
print("132 subspace frame?", any(L(list(U)) == f132 for U in [()]+[(a,) for a in range(1,8)]+list(itertools.combinations(range(1,8),2))+[tuple(range(1,8))]),
      " excess of 132:", d(L0, f132) + d(f132, LF) - m)
chk = dict(a0_in_XS=f57 in XS, phia0_in_XT=f134 in XT, astar_in_YS=f71 in YS, phiastar_in_YT=f47 in YT,
           mu_in_YS=f55 in YS, phimu_in_XT=f51 in XT, nu_in_XS=f79 in XS, phinu_in_YT=f71 in YT, cost=d(f57, f71))
print(chk)
# find a symplectic phi with phi(57)=134, phi(71)=47, phi(55)=51, phi(79)=71 ; DFS on images of basis
want = [(f57, f134), (f71, f47), (f55, f51), (f79, f71)]
def apply(cols, v):
    r = 0
    for b in range(6):
        if (v >> b) & 1: r ^= cols[b]
    return r
sols = []
zimgs = [v for v in f47 if v]
for zc in itertools.permutations(zimgs, 3):
    if len(span(zc)) != 8: continue
    # x images: om(x_i, z_j) = delta_ij, om(x_i, x_j) = 0
    cand = [[v for v in range(64) if all(om(v, zc[j]) == (i == j) for j in range(3))] for i in range(3)]
    for xc in itertools.product(*cand):
        if om(xc[0], xc[1]) or om(xc[0], xc[2]) or om(xc[1], xc[2]): continue
        cols = list(xc) + list(zc)   # bit b<3 -> x_b, bit 3+j -> z_j
        if all(frozenset(apply(cols, v) for v in A) == B for A, B in want):
            sols.append(cols)
print("symplectic phi satisfying all four frame conditions:", len(sols))
cols = sols[0]
fmt = lambda v: "(%s|%s)" % ("".join(str((v >> i) & 1) for i in range(3)), "".join(str((v >> (3 + i)) & 1) for i in range(3)))
print("one phi: images of (100|000),(010|000),(001|000),(000|100),(000|010),(000|001):", [fmt(c) for c in cols])
print("phi(134) =", names.get(frozenset(apply(cols, v) for v in f134), "?"), "in every solution:",
      all(frozenset(apply(c, v) for v in f134) == f132 for c in sols))
# interval structure: midpoints
print("midpoints of 57..71:", sorted(names.get(Z, "?") for Z in allL if d(f57, Z) == 1 and d(Z, f71) == 1))
# floor check: p_T = (u_T|u_T) (frames.py convention); is (u_T|0) in X_T frames?
pT0 = uT; pT1 = uT | (uT << m); qS = uS << m
print("(u_T|0) in every X_T frame:", all(pT0 in Z for Z in XT), "; in some:", any(pT0 in Z for Z in XT),
      "| (u_T|u_T) in every X_T frame:", all(pT1 in Z for Z in XT), "| (0|u_S) in every Y_S frame:", all(qS in Z for Z in YS))
ph = lambda v: apply(cols, v)
print("twisted pairing omega(p_T, phi q_S) =", om(pT1, ph(qS)), "; omega(p_S, phi^-1 q_T) via omega(phi p_S, q_T) =", om(ph(uS | (uS << m)), uT << m))
