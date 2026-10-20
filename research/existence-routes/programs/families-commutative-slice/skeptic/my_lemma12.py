# Lemma 12' and the Corollary 15 hyperplane claim, brute force with own code.
import geo, sys
from itertools import combinations
for m in range(4, int(sys.argv[1]) + 1):
    F = frozenset(z for z in range(1, 1 << m) if (z & 1) and (z & 2))
    assert geo.is_form_support(F, m) and geo.trace(F, m) == 0
    # odd (m-2)-flats avoiding 0 with exactly one basis vector (Lemma 12' asks for this), and all odd ones
    ev, od = geo.even_odd_flats(m, m - 2)
    od1 = [G for G in od if sum(1 for j in range(m) if (1 << j) in G) == 1]
    od = [frozenset(G) for G in od]; od1 = [frozenset(G) for G in od1]
    tot = 0; unc = 0; unc1 = 0
    for d in range(0, m - 1):
        subs = [A for A in geo.flats_avoiding_0(m, d) if set(A) <= F]
        for A in subs:
            A = frozenset(A)
            c_any = any(A <= G for G in od); c_one = any(A <= G for G in od1)
            if d <= m - 3:
                tot += 1; unc += (not c_any); unc1 += (not c_one)
            else:
                print("  m=%d: F itself (d=%d) inside an odd (m-2)-flat: %s (must be False)" % (m, d, c_any))
    print("m=%d: subflats of F of dim <= m-3: %d; not inside any odd (m-2)-flat avoiding 0: %d; not inside one with exactly one basis vector: %d"
          % (m, tot, unc, unc1))
