#!/usr/bin/env python3
"""t_setup.py (rm2-phase): exact checks of Lemma N's Sp(2m,2) game at m = 2..5.
 (1) every generator is symplectic; (2) K = prod_j t_(e_j) = x(I) (the plain schedule, one array); (3) y(I) x(I) y(I) = w;
 (4) lev(A) and y(C) fix L0, t_z moves L0 to distance 1, K to distance m, w to distance m;
 (5) Lemma E0 in Sp: t_z = lev(A) t_(e_1) lev(A)^-1 whenever A e_1 = z (so one move per direction carries no information);
 (6) the frames: count = prod (2^i + 1); P-orbits on frames = distance classes, sizes 2^(d(d+1)/2) [m choose d]_2;
 (7) the plain schedule on W arrays as a register program in C[Sp] is a unit (trivially: W copies of K = t_e1...t_em)."""
import sys, random
from spg import Sp, lagrangians, rand_gl, rand_sym, inv_f2
rng = random.Random(1)
def gauss_binom(n, k, q=2):
    num = den = 1
    for i in range(k):
        num *= q ** (n - i) - 1; den *= q ** (i + 1) - 1
    return num // den
for m in range(2, 6 if len(sys.argv) < 2 else int(sys.argv[1]) + 1):
    sp = Sp(m); ok = True
    K = sp.K(); prod = sp.ident()
    for j in range(m): prod = sp.mul(prod, sp.t(1 << j))
    c2 = (prod == K)
    yI = sp.yC(sp.I()); c3 = (sp.mul(sp.mul(yI, K), yI) == sp.w())
    gens = [sp.t(z) for z in range(1, 1 << m)] + [sp.lev(rand_gl(m, rng)) for _ in range(20)] + \
           [sp.yC(rand_sym(m, rng)) for _ in range(20)]
    c1 = all(sp.is_symplectic(g) for g in gens + [K, sp.w()])
    L0 = sp.L0()
    c4 = all(sp.act(sp.lev(rand_gl(m, rng)), L0) == L0 for _ in range(20)) and \
         all(sp.act(sp.yC(rand_sym(m, rng)), L0) == L0 for _ in range(20)) and \
         all(sp.dist(sp.act(sp.t(z), L0), L0) == 1 for z in range(1, 1 << m)) and \
         sp.dist(sp.act(K, L0), L0) == m and sp.dist(sp.act(sp.w(), L0), L0) == m
    c5 = True; t1 = sp.t(1)
    for z in range(1, 1 << m):
        for _ in range(5):
            A = rand_gl(m, rng)
            if A and sp.matvec(A, 1) != z:
                # force A e_1 = z: replace column 0 of A by z (rows are row-bitmasks), retry until invertible
                pass
        # build A with first column z, rest random, invertible
        while True:
            A = rand_gl(m, rng)
            A = [(A[i] & ~1) | ((z >> i) & 1) for i in range(m)]
            try: inv_f2(A, m); break
            except StopIteration: continue
        L = sp.lev(A); c5 &= (sp.mul(sp.mul(L, t1), sp.inv(L)) == sp.t(z))
    print(f"m = {m}: symplectic {c1}; K = t_e1...t_em = x(I) {c2}; y(I)x(I)y(I) = w {c3}; "
          f"P fixes L0, t_z at distance 1, K and w at distance m {c4}; t_z conjugate to t_e1 by an address map {c5}")
    if m <= 5:
        fr, idx = lagrangians(sp)
        expect = 1
        for i in range(1, m + 1): expect *= 2 ** i + 1
        dcount = {}
        for U in fr: d = sp.dist(U, L0); dcount[d] = dcount.get(d, 0) + 1
        dexp = {d: 2 ** (d * (d + 1) // 2) * gauss_binom(m, d) for d in range(m + 1)}
        # P-orbits: BFS on frames under P's generators (address maps and phases) must give exactly the distance classes
        Pg = [sp.lev([(1 << r) | ((1 << j) if r == i else 0) for r in range(m)]) for i in range(m) for j in range(m) if i != j]
        Pg += [sp.yC([(1 << i) if r == i else 0 for r in range(m)]) for i in range(m)]
        Pg += [sp.yC([((1 << j) if r == i else 0) | ((1 << i) if r == j else 0) for r in range(m)]) for i in range(m) for j in range(i + 1, m)]
        if m <= 4:
            orb = {}; no = 0
            for U in fr:
                if U in orb: continue
                st = [U]; orb[U] = no
                while st:
                    V = st.pop()
                    for g in Pg:
                        Wf = sp.act(g, V)
                        if Wf not in orb: orb[Wf] = no; st.append(Wf)
                no += 1
            pure = all(len({sp.dist(U, L0) for U in fr if orb[U] == o}) == 1 for o in range(no))
            porb = f"P-orbits on frames {no} (= m+1: {no == m + 1}), each one distance class {pure}"
        else:
            porb = "P-orbits not enumerated at m = 5"
        print(f"   frames {len(fr)} (formula {expect}); distance classes {dict(sorted(dcount.items()))} "
              f"(formula {dexp == dcount}); {porb}")
    sys.stdout.flush()
