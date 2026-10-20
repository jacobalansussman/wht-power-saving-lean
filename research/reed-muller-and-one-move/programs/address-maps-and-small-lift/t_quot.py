"""t_quot.py (rm-theory, claim Q1): single-array cost in the Gold quotient games of the trace placement.
Quotient of C[Sym_m x| Z_n] by the Singer-invariant part of Alt_m outside one cyclotomic class 1 + 2^i: the group
F_q (+) F_q x| Z_n, moves [(z, z^d)] (d = 1 + 2^i), target [(1, 0)]. A single array served by r moves needs
y_1..y_r in F_q^* with sum y = 1 and sum y^d = 0 (products of moves in an abelian normal subgroup; free conjugation
multiplies y by alpha^k). Reports the least r (2, 3, ...) at which a solution exists.
Two quotient classes at once (d1, d2): sum y = 1, sum y^d1 = sum y^d2 = 0."""
import itertools, sys
from gf import GF

def least_r(F, ds, rmax=6):
    q = F.q; elts = list(range(1, q))
    for r in range(1, rmax + 1):
        # meet in the middle for r >= 4 would be faster; small m only
        if r <= 3:
            for ys in itertools.combinations_with_replacement(elts, r):
                s = 0
                for y in ys: s ^= y
                if s != 1: continue
                if all(_pw_sum(F, ys, d) == 0 for d in ds): return r, ys
        else:
            # r = 4: pairs table
            from collections import defaultdict
            tab = defaultdict(list)
            for a, b in itertools.combinations_with_replacement(elts, 2):
                key = (a ^ b,) + tuple(F.pw(a, d) ^ F.pw(b, d) for d in ds)
                tab[key].append((a, b))
            for key, lst in tab.items():
                want = (key[0] ^ 1,) + tuple(key[1:])
                if want in tab:
                    return 4, lst[0] + tab[want][0]
            return None
    return None

def _pw_sum(F, ys, d):
    s = 0
    for y in ys: s ^= F.pw(y, d)
    return s

for m in [int(a) for a in sys.argv[1:]] or [4, 5, 6, 7]:
    F = GF(m)
    for i in range(1, m // 2 + 1):
        d = 1 + (1 << i)
        print(f"m={m} class d={d}: least r = {least_r(F, [d])}")
    if m >= 5:
        ds = [1 + (1 << i) for i in range(1, m // 2 + 1)]
        print(f"m={m} all classes {ds} (the full Sym_m quotient = the slice target I): least r = {least_r(F, ds, 4)}")

# two-class quotients (proper when m >= 6)
if len(sys.argv) > 1 and sys.argv[1] == "two":
    pass
