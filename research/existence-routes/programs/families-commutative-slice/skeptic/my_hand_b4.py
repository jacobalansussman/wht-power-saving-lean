# Independent exact check of the hand proof of Theorem 17 (weights 18/17 on the 7 even planes of E*).
from fractions import Fraction as Fr
import geo
m = 4
ev, od = geo.even_odd_flats(m, 2)
od = [frozenset(G) for G in od]
assert len(od) == 56 and len(ev) == 49
Estar = [z for z in range(1, 16) if bin(z).count("1") % 2 == 0]
assert len(Estar) == 7
planes = [frozenset(P) for P in ev if set(P) <= set(Estar)]
print("even planes inside E*:", [sorted(P) for P in planes])
assert len(planes) == 7
for P in planes:
    assert geo.is_form_support(P, m) and geo.trace(P, m) == 0
def avg_slack(U):
    Gs = [G for G in od if U <= G]
    assert Gs
    coef = {z: Fr(0) for z in range(1, 16)}
    for G in Gs:
        for z in G:
            coef[z] += Fr(1, len(Gs))
    return coef, Fr(-1)   # average of (sum_G c - 1)
total = {z: Fr(0) for z in range(1, 16)}; const = Fr(0); wsum = 0
for P in planes:
    w = 18 if 15 in P else 17
    wsum += w
    # RHS_P - LHS_P
    for x in P:
        total[x] -= w
    for x in P:
        for y in P:
            co, k = avg_slack(frozenset({x, y}))
            for z in co:
                total[z] += w * co[z]
            const += w * k
print("weights sum", wsum, "types A", sum(1 for P in planes if 15 in P))
print("coefficients of c_z in sum_P w_P (RHS_P - LHS_P):", sorted(set(total.values())))
print("constant:", const)
print("=> 0 <= %s * sum c + (%s)  =>  sum c >= %s" % (total[1], const, -const / total[1]) if len(set(total.values())) == 1 else "NOT uniform")
