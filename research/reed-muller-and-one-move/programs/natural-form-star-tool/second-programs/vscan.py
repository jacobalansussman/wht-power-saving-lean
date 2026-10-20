"""vscan.py (rm2-exact, independent of xbase / gen_tables / scanrc): verify every line of scanrc output.

Lines: "S k1,..,kW X F phi psi bound sF G_H G_H' G_psi0 G_x0 G_x1 G_psi1"  (ten-piece exclusion)
       "S k1,..,kW Y F phi psi c1 c2 2|Psi| sum G1 G2 G3"                (star-lemma exclusion)
       "S k1,..,kW D"                                                    (deletion fails)
       "S k1,..,kW P"                                                    (passes; nothing to verify)
Model (block 5, one move per nonzero direction, array l owed T_l = mult(alpha^(2 k_l)), alpha a root of x^5+x^2+1):
  - own GF(32) arithmetic (shift-and-reduce), own trace; the coordinate basis is checked to be self-dual;
  - a form index q is decoded by definition: Q(z) = sum_{p=(j<=k)} q_p z_j z_k (pairs in lexicographic order), and
    its support is mapped to field elements x = sum_j z_j f_j;
  - parity of Q at T_k: sum_j Q(alpha^k f_j) mod 2 (T_k = g g^T with g = mult(alpha^k), symmetric in this basis);
  - slack s_Q = |supp Q| - sum_l parity(Q, T_{k_l})   (rank E[D_Q, D_Q], Theorem 10 with a target per array);
  - X: the ten pieces are rebuilt from F, phi, psi; they must partition supp F x supp F; each piece R x C gets
    min(|R|, |C|, s_G) with R u C inside supp G (G as named); the sum must be < s_F.
  - D: some form must have s_Q < 0 (searched over all forms).
usage: vscan.py <basis as 5 ints> <scan output files...>
"""
import sys
import numpy as np

POLY = 0b100101
def gmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & 32: a ^= POLY
    return r
def tr(a):
    s, y = 0, a
    for _ in range(5): s ^= y; y = gmul(y, y)
    return s

basis = [int(x) for x in sys.argv[1].split(",")]
for i in range(5):
    for j in range(5):
        assert tr(gmul(basis[i], basis[j])) == (i == j)
for beta in range(1, 32):
    Mt = [[tr(gmul(beta, gmul(basis[j], basis[k]))) for k in range(5)] for j in range(5)]
    assert all(Mt[j][k] == Mt[k][j] for j in range(5) for k in range(5))
pw = [1]
for _ in range(30): pw.append(gmul(pw[-1], 2))
assert gmul(pw[-1], 2) == 1 and len(set(pw)) == 31
coord = {}
for z in range(1, 32):
    x = 0
    for j in range(5):
        if z >> j & 1: x ^= basis[j]
    assert all(tr(gmul(x, basis[j])) == (z >> j & 1) for j in range(5))
    coord[x] = z
PAIRS = [(j, k) for j in range(5) for k in range(j, 5)]
NF = 1 << 15
supp = []; par = np.zeros((NF, 31), dtype=np.int64); size = np.zeros(NF, dtype=np.int64)
for q in range(NF):
    S = set()
    for x, z in coord.items():
        v = 0
        for p, (j, k) in enumerate(PAIRS):
            if q >> p & 1: v ^= (z >> j & 1) & (z >> k & 1)
        if v: S.add(x)
    supp.append(frozenset(S)); size[q] = len(S)
    for k in range(31):
        par[q, k] = sum(1 for f in basis if gmul(pw[k], f) in S) % 2

def lin(chi, x):
    z = coord[x]; return bin(chi & z).count("1") & 1

nX = nD = nP = nY = 0; bad = 0
for fn in sys.argv[2:]:
    for line in open(fn):
        p = line.split()
        if not p or p[0] != "S": continue
        ks = [int(x) for x in p[1].split(",")]
        cnt = np.zeros(31, dtype=np.int64)
        for k in ks: cnt[k % 31] += 1
        r = par @ cnt; s = size - r
        if p[2] == "P": nP += 1; continue
        if p[2] == "D":
            if not (s[1:] < 0).any(): bad += 1; print("BAD D", line.strip())
            nD += 1; continue
        if p[2] == "Y":
            # star lemma, k = 0 (e3-slice6 7.2): E[D_F, D_F] nonsingular (s_F = |D_F|); Psi = D_F n {phi=c1, psi=c2};
            # Q_j = the rest of the three affine hyperplanes through that 3-flat; 2|Psi| <= sum_j rank E[Psi u Q_j].
            F, ph, ps, c1, c2, lhs, rhs = map(int, p[3:10]); Gs = list(map(int, p[10:13]))
            SF = supp[F]
            assert s[F] == len(SF), "E[D_F,D_F] not known to be nonsingular"
            hyp = [lambda x: lin(ph, x) == c1, lambda x: lin(ps, x) == c2, lambda x: lin(ph ^ ps, x) == (c1 ^ c2)]
            Psi = {x for x in SF if lin(ph, x) == c1 and lin(ps, x) == c2}
            parts = [{x for x in SF if hp(x)} - Psi for hp in hyp]
            assert Psi.isdisjoint(parts[0] | parts[1] | parts[2])
            assert sum(len(q) for q in parts) + len(Psi) == len(SF) and (Psi | parts[0] | parts[1] | parts[2]) == SF
            tot = 0
            for q, G in zip(parts, Gs):
                U = Psi | q; b = len(U)
                if G >= 0 and G != F and U <= supp[G]: b = min(b, int(s[G]))
                tot += b
            if not (2 * len(Psi) > tot and tot <= rhs): bad += 1; print("BAD Y", line.strip(), 2 * len(Psi), tot)
            nY += 1; continue
        assert p[2] == "X"
        F, ph, ps, bnd, sF = map(int, p[3:8]); Gs = list(map(int, p[8:14]))
        SF = supp[F]
        assert s[F] == sF, ("sF mismatch", line)
        H0 = {x for x in SF if lin(ph, x) == 0}; H1 = SF - H0
        K = {x for x in H0 if lin(ps, x) == 0}; L = H0 - K
        K1 = {x for x in H1 if lin(ps, x) == 0}; L1 = H1 - K1
        pieces = [(H0, H0, Gs[0]), (H1, H1, Gs[1]),
                  (K, K1, Gs[2]), (K1, K, Gs[2]), (K, L1, Gs[3]), (L1, K, Gs[3]),
                  (L, K1, Gs[4]), (K1, L, Gs[4]), (L, L1, Gs[5]), (L1, L, Gs[5])]
        cells = set(); total = 0
        for (R, C, G) in pieces:
            for x in R:
                for y in C:
                    assert (x, y) not in cells; cells.add((x, y))
            b = min(len(R), len(C))
            if R and C and G >= 0 and G != F and (R | C) <= supp[G]:
                b = min(b, int(s[G]))
            total += b
        ok = len(cells) == len(SF) ** 2 and total < sF and total <= bnd
        if not ok: bad += 1; print("BAD X", line.strip(), total, sF)
        nX += 1
print("verified: X %d, Y %d, D %d, (P %d not checked), BAD %d" % (nX, nY, nD, nP, bad))
