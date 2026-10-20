"""struct_check.py (rm-singer): exact checks of the canonical Singer formulation at m = 2..5.
 (1) sigma is symmetric in self-dual coordinates; so the virtual targets sigma^-k sigma^-kT = sigma^-2k are the
     trace forms Tr(c x y), c = alpha^-2k.
 (2) LINEAR-PART THEOREM: for every binary quadratic form Q and every channel k,
       <Q, sigma^-k sigma^-kT> = Tr(a_Q alpha^-k),  a_Q = linear coefficient of Q's univariate polynomial;
     and with move direction z0 = 1 the move is odd on channel k iff Q(alpha^-k) = 1.
 (3) the vectors (Q(alpha^-k))_k form a cyclic code of dimension m(m+1)/2 (punctured RM(2,m) without constant).
 (4) Z_n-orbits of forms: sizes, stabilisers; for each orbit |K|, |J|; odd orbits <-> purely quadratic q.
 (5) the plain identity sum_i Q(d w_i) = Tr(a_Q d) for all d (the plain schedule in Singer form).
"""
import sys
from collections import Counter
from sfield import Singer

for m in [int(a) for a in sys.argv[1:]] or [2, 3, 4, 5]:
    S = Singer(m); F = S.F; n = S.n
    sym = S.sigma == S.tp(S.sigma)
    print(f"m={m} n={n} self-dual basis {S.w} ; sigma symmetric: {sym}")
    assert sym
    # (1) T_k is the matrix of Tr(alpha^-2k x y)
    for k in range(n):
        c = F.a_exp(-2 * k)
        T = [[F.tr(F.mul(c, F.mul(S.w[i], S.w[j]))) for j in range(m)] for i in range(m)]
        assert T == S.Tk(k)
    print("  (1) virtual targets = trace forms Tr(alpha^-2k xy): ok")
    # (2), (5)
    bad = 0; nforms = 0; orbit_of = {}; orbits = []
    code = set()
    for Q in S.forms():
        nforms += 1
        a = S.linpart(Q)
        K, J = S.KJ(Q)
        Jlin = frozenset(k for k in range(n) if F.tr(F.mul(a, F.a_exp(-k))))
        if J != Jlin: bad += 1
        for d in range(1, F.q):
            if (sum(S.Qval(Q, F.mul(d, wi)) for wi in S.w) & 1) != F.tr(F.mul(a, d)): bad += 1
        code.add(tuple(int(k in K) for k in range(n)))
    print(f"  (2,5) linear-part theorem over {nforms} forms x {n} channels and plain identity: violations {bad}")
    assert bad == 0
    # (3) cyclic code
    cw = sorted(code)
    rot = lambda v: v[-1:] + v[:-1]
    cyc = all(rot(v) in code for v in cw)
    import math
    print(f"  (3) code size {len(code)} = 2^{math.log2(len(code)):.0f} (expect 2^{S.N}); cyclic: {cyc}")
    assert len(code) == 2 ** S.N and cyc
    # (4) orbits under Q -> Q o sigma^{-1}: realised on the K-vector by rotation (Q determined by its values on F*)
    seen = set(); stats = Counter()
    for v in cw:
        if v in seen: continue
        orb = [v]; w = rot(v)
        while w != v: orb.append(w); w = rot(w)
        seen |= set(orb)
        # find the form: recompute J via linear part from a representative (rebuild Q from values is not needed:
        # J depends only on the linear part; compute a from values directly: a = sum_{y} Q(y) y^-1)
        a = 0
        for k in range(n):
            if v[k]: a ^= F.a_exp(k)   # y = alpha^-k, y^-1 = alpha^k
        J = sum(F.tr(F.mul(a, F.a_exp(-k))) for k in range(n))
        stats[(len(orb), sum(v), J)] += 1
    print("  (4) orbits (size, |K|, |J|): count")
    for key in sorted(stats):
        print("     ", key, stats[key])
    odd = sum(c for (sz, K, J), c in stats.items() if J > 0)
    print(f"      odd orbits: {odd} (expect 2^(N-m) = {2 ** (S.N - m)}: one per purely quadratic q)")
print("ALL OK")
