"""vcert.py (rm2-exact, independent of xbase/rcdp): verify a rank-cover certificate for the multi-target system with
one move per nonzero direction and canonical Singer targets T_k = mult(alpha^(2k)), alpha a root of x^5 + x^2 + 1.

Everything is recomputed in FIELD terms from scratch:
  - GF(32) by shift-and-reduce multiplication (no log tables);
  - the coordinate basis f_1..f_5 (taken from the certificate header) is checked to be self-dual for Tr(xy);
  - a direction with coordinate vector z (bit j = coordinate j) is the field element x = sum_j z_j f_j (checked:
    Tr(x f_j) = z_j);
  - a form is given only as its support S (a set of nonzero field elements); it is checked to be a quadratic
    function vanishing at 0 (all third derivatives vanish);
  - its parity at target mult(beta^2) is sum_j [beta f_j in S] mod 2 (T = g g^T with g = mult(beta), which is symmetric
    in a self-dual basis -- checked);
  - slack s = |S| - #odd targets;  piece bound checked against min(|R|,|C|) or s_G with R u C inside S_G;
  - the pieces must partition S_F x S_F and their bounds must sum to less than s_F.
usage: vcert.py <certfile> <k1,k2,...>
"""
import sys
POLY = 0b100101
def gmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & 32: a ^= POLY
    return r
def gpow(a, e):
    r = 1
    for _ in range(e): r = gmul(r, a)
    return r
def tr(a):
    s, y = 0, a
    for _ in range(5): s ^= y; y = gmul(y, y)
    assert s in (0, 1); return s
ALPHA = 2
assert gpow(ALPHA, 31) == 1 and all(gpow(ALPHA, e) != 1 for e in range(1, 31))

def main():
    lines = open(sys.argv[1]).read().split("\n")
    ks = [int(x) for x in sys.argv[2].split(",")]
    hdr = [L for L in lines if L.startswith("# case")][0]
    basis = eval(hdr.split("basis")[1].split("poly")[0])
    assert len(basis) == 5
    for i in range(5):
        for j in range(5):
            assert tr(gmul(basis[i], basis[j])) == (1 if i == j else 0), "basis not self-dual"
    def elem(z):
        x = 0
        for j in range(5):
            if z >> j & 1: x ^= basis[j]
        for j in range(5): assert tr(gmul(x, basis[j])) == (z >> j & 1)
        return x
    def mask_to_set(mask):
        return frozenset(elem(i + 1) for i in range(31) if mask >> i & 1)
    # mult(beta) symmetric in this basis
    for beta in range(1, 32):
        Mt = [[tr(gmul(beta, gmul(basis[j], basis[k]))) for k in range(5)] for j in range(5)]
        assert all(Mt[j][k] == Mt[k][j] for j in range(5) for k in range(5))
    betas = [gpow(ALPHA, k % 31) for k in ks]
    def is_quadratic(S):
        Q = lambda x: 1 if x in S else 0
        if 0 in S: return False
        for x in range(32):
            for y in range(x, 32):
                for z in range(y, 32):
                    v = Q(x ^ y ^ z) ^ Q(x ^ y) ^ Q(x ^ z) ^ Q(y ^ z) ^ Q(x) ^ Q(y) ^ Q(z) ^ Q(0)
                    if v: return False
        return True
    def slack(S):
        assert is_quadratic(S), "support is not a quadratic form"
        r = 0
        for b in betas:
            r += sum(1 for f in basis if gmul(b, f) in S) % 2
        return len(S) - r
    PAIRS = [(j, k) for j in range(5) for k in range(j, 5)]
    def form_support(q):
        S = set()
        for z in range(1, 32):
            v = 0
            for p, (j, k) in enumerate(PAIRS):
                if q >> p & 1: v ^= (z >> j & 1) & (z >> k & 1)
            if v: S.add(elem(z))
        return frozenset(S)
    i = 0; nver = 0
    while i < len(lines):
        if lines[i].startswith("CERT"):
            F = int(lines[i].split()[1]); pieces = []; i += 1
            while lines[i] != "END":
                p = lines[i].split(); pieces.append((int(p[1]), int(p[2]), int(p[3]), int(p[4]), int(p[5]))); i += 1
            SF = form_support(F); sF = slack(SF)
            cells = set(); total = 0
            for (R, C, b, how, G) in pieces:
                Rs, Cs = mask_to_set(R), mask_to_set(C)
                assert Rs <= SF and Cs <= SF and Rs and Cs
                for x in Rs:
                    for y in Cs:
                        assert (x, y) not in cells, "overlap"; cells.add((x, y))
                if how == 0:
                    assert b >= min(len(Rs), len(Cs))
                else:
                    SG = form_support(G)
                    assert (Rs | Cs) <= SG, "piece not inside S_G"
                    assert b >= slack(SG), "bound below slack"
                total += b
            assert len(cells) == len(SF) ** 2, "not a partition"
            print("VERIFIED form %d: |S_F| = %d, s_F = %d, %d pieces, bound sum %d < s_F: %s" % (F, len(SF), sF, len(pieces), total, total < sF))
            assert total < sF
            nver += 1
        i += 1
    print("certificates verified:", nver)

if __name__ == "__main__":
    main()
