# Near-uniform multisets at block 6 (e3-slice6): U_6 plus extra moves Z (a multiset of directions), W arrays.
# Theorem C applies if some even 4-flat F (avoiding 0, even number of basis vectors) avoids Z and every affine
# hyperplane A of F lies in an odd 4-flat G (avoiding 0, odd number of basis vectors) with slack
# 16 + #(Z in G) - W <= 3.   (E[F,F] is then a 16 x 16 nonsingular matrix with all hyperplane ranks <= 3.)
import itertools, sys
m = 6; N = 1 << m
basis = {1 << j for j in range(m)}
def span(vs):
    S = {0}
    for v in vs: S |= {x ^ v for x in S}
    return S
subs4 = set()
for vs in itertools.combinations(range(1, N), 4):
    S = frozenset(span(vs))
    if len(S) == 16: subs4.add(S)
flats4 = set()
for S in subs4:
    for t in range(N):
        Fl = frozenset(x ^ t for x in S)
        if 0 not in Fl: flats4.add(Fl)
flats4 = list(flats4)
par = {Fl: len(Fl & basis) % 2 for Fl in flats4}
even4 = [Fl for Fl in flats4 if par[Fl] == 0]; odd4 = [Fl for Fl in flats4 if par[Fl] == 1]
def hyperplanes(Fl):
    pts = sorted(Fl); a = pts[0]; D = [p ^ a for p in pts]
    out = []
    for h in range(1, N):
        for b in (0, 1):
            A = frozenset(p for p, d in zip(pts, D) if bin(h & d).count("1") % 2 == b)
            if len(A) == 8: out.append(A)
    return list(set(out))
H_of = {Fl: hyperplanes(Fl) for Fl in even4}
odd_through = {}
for Fl in even4:
    for A in H_of[Fl]:
        if A not in odd_through: odd_through[A] = [G for G in odd4 if A <= G]
print("4-flats avoiding 0:", len(flats4), "even", len(even4), "odd", len(odd4))
good_F = [Fl for Fl in even4 if all(odd_through[A] for A in H_of[Fl])]
print("even 4-flats all of whose hyperplanes lie in an odd 4-flat avoiding 0:", len(good_F))
def excluded(Z, W):
    Zs = list(Z)
    for Fl in good_F:
        if any(z in Fl for z in Zs): continue
        if all(any(16 + sum(1 for z in Zs if z in G) - W <= 3 for G in odd_through[A]) for A in H_of[Fl]):
            return True
    return False
dirs = range(1, N)
for W, tmax in ((13, 2), (14, 3), (15, 3), (16, 3)):
    for t in range(0, tmax + 1):
        tot = exc = 0; bad = []
        for Z in itertools.combinations_with_replacement(dirs, t):
            tot += 1
            if excluded(Z, W): exc += 1
            elif len(bad) < 3: bad.append(Z)
        M = 63 + t
        print(f"W={W} t={t} (M={M}, moves/array {M/W:.3f}, above bar: {M/W <= 6**0.9}): {exc}/{tot} multisets U_6+Z excluded; not covered e.g. {bad}", flush=True)
