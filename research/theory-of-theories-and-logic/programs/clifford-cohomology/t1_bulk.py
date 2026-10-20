# t1_bulk.py -- meta-cohomology test 1: the bulk (Clifford QCA) level, exact.
# Tests, on random Clifford QCAs built from translation-invariant layers (m = 1, 2, 3):
#  (a) gamma(v, w) = sum_{d > d'} <M_d v, M_d' w> is identically 0 (prediction: the lag identities kill it);
#  (b) the naive ordered collapse Pi_0 is a valid one-column Clifford (images Hermitian, canonical commutation);
#  (c) Pi_0 is multiplicative: Pi_0(k l) == Pi_0(k) Pi_0(l), exactly with signs (the 2-cocycle question);
#  (d) the odd-ring restriction Psi_L (L = 7, 9, 11, 13, 15) is a valid Clifford, independent of L, and multiplicative;
#  (e) Psi_L versus Pi_0: they differ by the character chi(u) = <M_odd u, M_even u> (predicted) -- record whether chi
#      is ever nonzero;
#  (f) Pi_0 is the identity on constants (one-column layers) and on translation-invariant Paulis;
#  (g) the PROOF 7.5 relation: E^-1 Dcz E Dcz^-1 = CZ chain on row 1 (bulk), and its collapse = Z_1.
import random, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cqca import *

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 20261014)


def onecol(m, gate_fn_of_bits):
    """layer with the same one-column gate list at every site; gate_fn_of_bits(j) -> list of gates."""
    return layer_qca(m, gate_fn_of_bits, reach=2)


def elementary(m):
    """list of (name, QCA) elementary letters for block size m."""
    L = []
    for q in range(m):
        L.append(("H%d" % q, layer_qca(m, lambda j, q=q: [g_H(bit(m, j, q))], 1)))
        L.append(("S%d" % q, layer_qca(m, lambda j, q=q: [g_S(bit(m, j, q))], 1)))
        L.append(("X%d" % q, layer_qca(m, lambda j, q=q: [g_X(bit(m, j, q))], 1)))
        L.append(("Z%d" % q, layer_qca(m, lambda j, q=q: [g_Z(bit(m, j, q))], 1)))
        L.append(("D%d" % q, layer_qca(m, lambda j, q=q: [g_rot(1 << bit(m, j, q), 0)], 1)))
        for d in (1, 2):
            L.append(("CZc%d%d+%d" % (q, q, d), layer_qca(m, lambda j, q=q, d=d: [g_CZ(bit(m, j, q), bit(m, j + d, q))], 4)))
        # row shift (cyclic rotation of one row in the finite game)
        rule = list(QCA.ident(m).rule)
        rule[q] = (1 << bit(m, 1, q), 0, 0)
        rule[m + q] = (0, 1 << bit(m, 1, q), 0)
        L.append(("sh%d" % q, QCA(m, rule)))
        rule = list(QCA.ident(m).rule)
        rule[q] = (1 << bit(m, -1, q), 0, 0)
        rule[m + q] = (0, 1 << bit(m, -1, q), 0)
        L.append(("shinv%d" % q, QCA(m, rule)))
        for p in range(m):
            if p == q:
                continue
            L.append(("CNOT%d%d" % (q, p), layer_qca(m, lambda j, q=q, p=p: [g_CNOT(bit(m, j, q), bit(m, j, p))], 1)))
            L.append(("CZ%d%d" % (q, p), layer_qca(m, lambda j, q=q, p=p: [g_CZ(bit(m, j, q), bit(m, j, p))], 1)))
            L.append(("E%d>%d+1" % (q, p), layer_qca(m, lambda j, q=q, p=p: [g_CNOT(bit(m, j, q), bit(m, j + 1, p))], 4)))
            L.append(("E%d>%d-1" % (q, p), layer_qca(m, lambda j, q=q, p=p: [g_CNOT(bit(m, j + 1, q), bit(m, j, p))], 4)))
            L.append(("CZc%d%d+1" % (q, p), layer_qca(m, lambda j, q=q, p=p: [g_CZ(bit(m, j, q), bit(m, j + 1, p))], 4)))
    return L


def inverse_by_power(Q, cap=2000):
    """not used for big words; inverse of a QCA by repeated multiplication is too slow in general."""
    raise NotImplementedError


def rand_word(Lett, n):
    return [random.choice(Lett) for _ in range(n)]


def word_qca(m, w):
    Q = QCA.ident(m)
    for name, q in w:
        Q = Q * q
    return Q


def check_letters_valid(m, Lett):
    bad = []
    for name, q in Lett:
        # validity of a QCA rule: images Hermitian; commutation with all shifts matches the canonical one
        r = q.radius()
        for i, P in enumerate(q.rule):
            if not is_herm(P):
                bad.append((name, "nonherm", i))
            for j, Qp in enumerate(q.rule):
                for s in range(-2 * r - 1, 2 * r + 2):
                    want = 1 if (s == 0 and abs(i - j) == m) else 0
                    if symp(P, shiftP(Qp, s, m)) != want:
                        bad.append((name, "comm", i, j, s))
    return bad


def conjugate(C):
    """complex conjugation twist c C c (c: i -> -i, X, Z fixed): images (x, z, k) -> (x, z, -k) after applying C to c(gen)."""
    return Cliff([(p[0], p[1], (-p[2]) % 4) for p in C.X], [(p[0], p[1], (-p[2]) % 4) for p in C.Z])


total = {"psi3_ne_conj_collapse": 0, "compared_psi1": 0, "gamma_nonzero": 0, "collapse_invalid": 0, "collapse_nonmult": 0, "psi_invalid": 0,
         "psi_L_dependent": 0, "psi_nonmult": 0, "psi_ne_collapse_chi": 0, "chi_nonzero": 0, "pairs": 0, "elems": 0}

for m in (1, 2, 3):
    Lett = elementary(m)
    bad = check_letters_valid(m, Lett)
    print("m=%d: %d elementary letters, invalid letter rules: %d" % (m, len(Lett), len(bad)))
    # (f) constants and Paulis
    for name, q in Lett:
        if q.radius() == 0:
            c = collapse(q)
            # one-column letter: collapse = the one-column Clifford itself (rule at site 0 read on one column)
            direct = Cliff([(site_factor(p, m, 0)[0], site_factor(p, m, 0)[1], p[2]) for p in q.rule[:m]],
                           [(site_factor(p, m, 0)[0], site_factor(p, m, 0)[1], p[2]) for p in q.rule[m:]])
            assert c == direct, name
    print("   (f) collapse = the letter itself on every one-column letter: PASS")
    nwords = {1: 60, 2: 60, 3: 40}[m]
    elems = []
    for t in range(nwords):
        w = rand_word(Lett, random.randint(1, 10))
        Q = word_qca(m, w)
        elems.append((w, Q))
    stats_chi = 0
    for w, Q in elems:
        total["elems"] += 1
        G = gamma_matrix(Q)
        if any(any(r) for r in G):
            total["gamma_nonzero"] += 1
        c = collapse(Q)
        if not c.valid():
            total["collapse_invalid"] += 1
        r = Q.radius()
        psis = []
        for L in (7, 9, 11, 13, 15):
            if L > 2 * r + 1:
                psis.append((L, ring_restrict(Q, L)))
        if not psis:
            continue
        for L, p in psis:
            if p is None or not p.valid():
                total["psi_invalid"] += 1
        # Psi_L depends only on L mod 4 (prediction: L = 1 mod 4 gives Pi_0, L = 3 mod 4 its complex conjugate)
        for res in (1, 3):
            grp = [p for L, p in psis if L % 4 == res]
            if any(p != grp[0] for p in grp):
                total["psi_L_dependent"] += 1
                print("   L-dependence within L = %d mod 4:" % res, [(L, p.key()) for L, p in psis])
        for L, p in psis:
            if L % 4 == 3 and p != conjugate(c):
                total["psi3_ne_conj_collapse"] += 1
        # (e) Psi = Pi_0 twisted by chi:  Psi(P_u) = (-1)^chi(u) Pi_0(P_u) on generators
        chi = chi_form(Q)
        if any(chi):
            total["chi_nonzero"] += 1
        p1 = [p for L, p in psis if L % 4 == 1]
        if not p1:
            continue
        p = p1[0]
        total["compared_psi1"] += 1
        tw = Cliff([(P[0], P[1], (P[2] + 2 * chi[i]) % 4) for i, P in enumerate(c.X)],
                   [(P[0], P[1], (P[2] + 2 * chi[c.n + i]) % 4) for i, P in enumerate(c.Z)])
        if p != tw:
            total["psi_ne_collapse_chi"] += 1
    # (c), (d): multiplicativity on random pairs
    for t in range(len(elems)):
        (w1, A), (w2, B) = random.choice(elems), random.choice(elems)
        AB = A * B
        total["pairs"] += 1
        if collapse(AB) != collapse(A) * collapse(B):
            total["collapse_nonmult"] += 1
        r = max(A.radius(), B.radius(), AB.radius())
        for L in (2 * r + 3, 2 * r + 5):          # one L of each class mod 4
            pa, pb, pab = ring_restrict(A, L), ring_restrict(B, L), ring_restrict(AB, L)
            if pab != pa * pb:
                total["psi_nonmult"] += 1
    print("   m=%d done: %d elements, %d pairs so far" % (m, total["elems"], total["pairs"]))

# (g) the relation of PROOF 7.5 at m = 2: E = row0 ^= row1 << 1 (x_{0,j} ^= x_{1,j-1}); Dcz = CZ(0,1) on every column
m = 2
E = layer_qca(m, lambda j: [g_CNOT(bit(m, j - 1, 1), bit(m, j, 0))], 4)
Ei = E  # CNOT layer is an involution
Dcz = layer_qca(m, lambda j: [g_CZ(bit(m, j, 0), bit(m, j, 1))], 1)
lhs = Ei * Dcz * E * Dcz
CZchain1 = layer_qca(m, lambda j: [g_CZ(bit(m, j, 1), bit(m, j + 1, 1))], 4)
print("(g) E^-1 Dcz E Dcz^-1 == CZ chain on row 1 (bulk):", lhs == CZchain1)
c = collapse(CZchain1)
Z1 = Cliff([(1, 0, 0), (2, 0, 2)], [(0, 1, 0), (0, 2, 0)])   # conjugation by Z on row 1: X_1 -> -X_1 (row 1 = bit 2)
print("    collapse(CZ chain row 1) == conjugation by Z_1:", c == Z1,
      "| collapse(E)^-1 collapse(Dcz) collapse(E) collapse(Dcz)^-1 == Z_1:",
      collapse(Ei) * collapse(Dcz) * collapse(E) * collapse(Dcz) == Z1)
print()
for k, v in total.items():
    print("%-22s %d" % (k, v))
ok = (total["gamma_nonzero"] == 0 and total["collapse_invalid"] == 0 and total["collapse_nonmult"] == 0
      and total["psi_invalid"] == 0 and total["psi_L_dependent"] == 0 and total["psi_nonmult"] == 0
      and total["psi_ne_collapse_chi"] == 0 and total["psi3_ne_conj_collapse"] == 0 and total["chi_nonzero"] == 0)
print("ALL PASS" if ok else "SOME CHECK FAILED")
