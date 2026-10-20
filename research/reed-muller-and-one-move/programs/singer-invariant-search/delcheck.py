"""delcheck.py (rm-singer): fibre data (|K|, |J|) per H-orbit for the two placements of the Singer group (sd: self-dual
basis, sigma symmetric; poly: polynomial basis), and the exact deletion test M|K| >= W|J| for given shapes."""
import sys, math
from collections import Counter
from sfield import Singer
for m in (3, 4, 5):
    for coords in ('sd', 'poly'):
        S = Singer(m, coords); F = S.F; n = S.n
        seen = set(); tab = Counter(); minratio = None
        for Q in S.forms():
            if not any(Q): continue
            v = tuple(S.Qval(Q, x) for x in range(1, F.q))
            if v in seen: continue
            orb = {tuple(S.Qval(Q, F.mul(F.a_exp(-e), x)) for x in range(1, F.q)) for e in range(n)}
            seen |= orb
            K = sum(S.Qval(Q, F.a_exp(-k)) for k in range(n))
            J = sum(S.pair(Q, S.Tk(k)) for k in range(n))
            tab[(len(orb), K, J)] += 1
            if J: 
                r = J / K
                minratio = r if minratio is None or r > minratio else minratio
        print(f"m={m} {coords}: sigma symmetric {S.sigma == S.tp(S.sigma)}; orbits (size,|K|,|J|):", dict(sorted(tab.items())))
        print(f"   deletion bound M/W >= max |J|/|K| = {minratio}")
