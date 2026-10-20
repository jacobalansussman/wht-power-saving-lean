"""x-transfer, 13 Oct 2026: frame geometry checks (standard library only; run with python3 -I).

Conventions as in x-dict/dictcheck.py and the hub's general-clifford-frames.tex:
a label vector is (x|z) in F_2^(2m), packed as the int (x << m) | z; omega((x,z),(x',z')) = x.z' + z.x'.
L_0 = {(0,z)} (start), L_F = {(z,z)} (end), L_U = {(u, u+v): u in U, v in U^perp} (subspace frames).
d(A,B) = m - dim(A cap B).

Checks:
 A. distance-one frames: per nonzero u exactly two, one subspace frame on a geodesic to L_F (d = m-1),
    one non-subspace frame at d(., L_F) = m.  m = 1..9 (1,022 at m = 9).
 B. excess e(Z) = d(L_0,Z) + d(Z,L_F) - m equals rank(beta_Z - dot|U) where U = X-projection of Z and beta_Z
    its Dehaene-De Moor form; table of frames by (r = d(L_0,Z), e) vs [m r]_2 * #sym(r, rank e).  m = 1..4.
 C. orbits on the 135 frames at m = 3: (i) free group = stabiliser of L_0 (Bruhat cells, by r);
    (ii) local Clifford group (H_j, S_j on each qubit): LC orbits; (iii) the problem symmetry group =
    stabiliser of both L_0 and L_F (the hub's K_G, = GL(m,2)); and the classification of frames opposite L_0
    ((x | Gamma x), "graph states with S twists") by congruence of Gamma + I.
"""
import itertools, sys

def popc(v): return bin(v).count("1")
def dot(a, b): return popc(a & b) & 1

def omega(m, p, q):
    mask = (1 << m) - 1
    return dot(p >> m, q & mask) ^ dot(p & mask, q >> m)

def span(vecs):
    s = {0}
    for v in vecs:
        if v not in s:
            s |= {w ^ v for w in s}
    return frozenset(s)

def dim_of(s): return len(s).bit_length() - 1

def L0(m): return span([1 << j for j in range(m)])                       # (0|e_j)
def LF(m): return span([((1 << j) << m) | (1 << j) for j in range(m)])     # (e_j|e_j)

def perp_in(m, U):  # U subset of F_2^m (set), its dot-perp
    return frozenset(v for v in range(1 << m) if all(dot(v, u) == 0 for u in U))

def LU(m, U):
    Up = perp_in(m, U)
    return frozenset((u << m) | (u ^ v) for u in U for v in Up)

def all_subspaces(m):
    seen = set()
    for k in range(m + 1):
        for basis in itertools.combinations(range(1, 1 << m), k):
            S = span(basis)
            if dim_of(S) == k:
                seen.add(S)
    return seen

def all_lagrangians(m):
    """All Lagrangians via (U, beta) (Dehaene-De Moor): Z = {(u, B u + v) : u in U, v in U^perp} with B a
    symmetric choice on a basis of U.  Deduplicated."""
    out = set()
    for U in all_subspaces(m):
        k = dim_of(U)
        basis = []
        for v in sorted(U):
            if v and v not in span(basis):
                basis.append(v)
        Up = perp_in(m, U)
        # images w_i of the basis vectors, chosen in F_2^m with w_i . u_j symmetric; enumerate all w tuples
        # modulo U^perp: pick representatives by trying all vectors and deduplicate the resulting sets.
        for ws in itertools.product(range(1 << m), repeat=k):
            ok = all(dot(ws[i], basis[j]) == dot(ws[j], basis[i]) for i in range(k) for j in range(i + 1, k))
            if not ok:
                continue
            gens = [(basis[i] << m) | ws[i] for i in range(k)] + [v for v in Up]
            Z = span(gens)
            if dim_of(Z) == m:
                out.add(Z)
    return out

def is_isotropic(m, Z):
    return all(omega(m, p, q) == 0 for p in Z for q in Z)

def dist(m, A, B): return m - dim_of(A & B)

def gauss_count(m, r):
    num = den = 1
    for i in range(r):
        num *= (1 << (m - i)) - 1
        den *= (1 << (i + 1)) - 1
    return num // den

def rank_f2(rows):
    rows = list(rows); r = 0
    while rows:
        p = rows.pop()
        if p == 0:
            continue
        r += 1
        lb = p & -p
        rows = [q ^ p if q & lb else q for q in rows]
    return r

def sym_rank_counts(r):
    """number of symmetric r x r matrices over F_2 by rank, brute force (r <= 4)."""
    cnt = {}
    pos = [(i, j) for i in range(r) for j in range(i, r)]
    for bits in range(1 << len(pos)):
        M = [0] * r
        for t, (i, j) in enumerate(pos):
            if bits >> t & 1:
                M[i] |= 1 << j; M[j] |= 1 << i
        k = rank_f2(M)
        cnt[k] = cnt.get(k, 0) + 1
    return cnt

def check_A(mmax):
    print("== A. distance-one frames from L_0, by distance to L_F")
    for m in range(1, mmax + 1):
        l0, lf = L0(m), LF(m)
        tally = {}
        for u in range(1, 1 << m):
            Hu = [z for z in range(1 << m) if dot(z, u) == 0]     # (0|z), z perp u
            w = next(z for z in range(1 << m) if dot(z, u) == 1)
            A = span(Hu + [u << m])                 # H_u + (u|0)
            B = span(Hu + [(u << m) | w])           # H_u + (u|w)
            sub = LU(m, span([u]))
            assert A != B and dist(m, l0, A) == 1 and dist(m, l0, B) == 1
            assert (A == sub) != (B == sub)
            for Z in (A, B):
                key = ("subspace" if Z == sub else "non-subspace", dist(m, Z, lf),
                       "u odd" if dot(u, u) else "u even")
                tally[key] = tally.get(key, 0) + 1
        tot = sum(tally.values())
        print(f"  m={m}: {tot} frames at distance 1 (2(2^m-1) = {2 * ((1 << m) - 1)});",
              "; ".join(f"{k[0]} at d(.,L_F)={k[1]} ({k[2]}): {v}" for k, v in sorted(tally.items())))

def beta_minus_dot_rank(m, Z):
    """Z Lagrangian.  U = X-projection; for u in U pick z with (u|z) in Z; beta(u,u') = z . u'.
    Return (dim U, rank of (beta - dot) on U)."""
    mask = (1 << m) - 1
    U = span([p >> m for p in Z])
    basis = []
    for v in sorted(U):
        if v and v not in span(basis):
            basis.append(v)
    zs = []
    for b in basis:
        zs.append(next(p & mask for p in Z if p >> m == b))
    k = len(basis)
    rows = []
    for i in range(k):
        row = 0
        for j in range(k):
            if dot(zs[i], basis[j]) ^ dot(basis[i], basis[j]):
                row |= 1 << j
        rows.append(row)
    # symmetric check
    for i in range(k):
        for j in range(k):
            assert (rows[i] >> j & 1) == (rows[j] >> i & 1)
    return k, rank_f2(rows)

def check_B(mmax):
    print("== B. excess e(Z) = d(L_0,Z)+d(Z,L_F)-m  versus  rank(beta_Z - dot) on the X-projection")
    symc = {r: sym_rank_counts(r) for r in range(mmax + 1)}
    for m in range(1, mmax + 1):
        Ls = all_lagrangians(m)
        l0, lf = L0(m), LF(m)
        subs = {LU(m, U) for U in all_subspaces(m)}
        bad = 0; table = {}
        for Z in Ls:
            assert is_isotropic(m, Z) and dim_of(Z) == m
            r, e = beta_minus_dot_rank(m, Z)
            exc = dist(m, l0, Z) + dist(m, Z, lf) - m
            if r != dist(m, l0, Z) or e != exc or ((e == 0) != (Z in subs)):
                bad += 1
            table[(r, e)] = table.get((r, e), 0) + 1
        formula = {(r, e): gauss_count(m, r) * c for r in range(m + 1) for e, c in symc[r].items()}
        print(f"  m={m}: {len(Ls)} Lagrangians; mismatches of r = d(L_0,Z), e = rank(beta-dot), e=0 <=> subspace: {bad};"
              f" table (r,e)->count equals [m r]_2 * #sym(r x r, rank e): {table == formula}")
        print("       " + ", ".join(f"(r={r},e={e}):{table[(r, e)]}" for (r, e) in sorted(table)))

def apply_sym(m, M, p):
    """M: list of 2m column images (ints) of basis vectors; basis index t is bit t of p."""
    out = 0
    for t in range(2 * m):
        if p >> t & 1:
            out ^= M[t]
    return out

def check_C():
    m = 3
    print("== C. orbits on the 135 Lagrangian frames at m = 3")
    Ls = sorted(all_lagrangians(m), key=lambda s: sorted(s))
    idx = {Z: i for i, Z in enumerate(Ls)}
    l0, lf = L0(m), LF(m)
    # generators as symplectic maps on packed vectors: bit t<m is z_t, bit m+t is x_t
    def H(j):      # swap x_j, z_j
        return lambda p: (p & ~((1 << j) | (1 << (m + j)))) | ((p >> j & 1) << (m + j)) | ((p >> (m + j) & 1) << j)
    def S(j):      # z_j += x_j
        return lambda p: p ^ (((p >> (m + j)) & 1) << j)
    def CX(c, t):  # x_t += x_c ; z_c += z_t
        return lambda p: p ^ (((p >> (m + c)) & 1) << (m + t)) ^ (((p >> t) & 1) << c)
    def CZ(a, b):  # z_a += x_b, z_b += x_a
        return lambda p: p ^ (((p >> (m + b)) & 1) << a) ^ (((p >> (m + a)) & 1) << b)
    def orbits(gens):
        seen = {}; orbs = []
        for Z in Ls:
            if Z in seen:
                continue
            o = {Z}; stack = [Z]
            while stack:
                Y = stack.pop()
                for g in gens:
                    W = frozenset(g(p) for p in Y)
                    if W not in o:
                        o.add(W); stack.append(W)
            for Y in o:
                seen[Y] = len(orbs)
            orbs.append(o)
        return orbs, seen
    for name, gens in [("free (CNOT, CZ, S; stabiliser of L_0)", [CX(a, b) for a in range(m) for b in range(m) if a != b] + [CZ(a, b) for a in range(m) for b in range(a + 1, m)] + [S(j) for j in range(m)]),
                       ("local Clifford (H_j, S_j)", [H(j) for j in range(m)] + [S(j) for j in range(m)])]:
        orbs, lab = orbits(gens)
        sizes = sorted(len(o) for o in orbs)
        print(f"  {name}: {len(orbs)} orbits, sizes {sizes}")
        if name.startswith("free"):
            assert all(len({dist(m, l0, Z) for Z in o}) == 1 for o in orbs)
            free_lab = lab
        else:
            lc_lab = lab
    # compare partitions
    pairs_same_free_diff_lc = sum(1 for A in Ls for B in Ls if free_lab[A] == free_lab[B] and lc_lab[A] != lc_lab[B])
    pairs_same_lc_diff_free = sum(1 for A in Ls for B in Ls if lc_lab[A] == lc_lab[B] and free_lab[A] != free_lab[B])
    print(f"  ordered pairs in one free orbit but two LC orbits: {pairs_same_free_diff_lc}; in one LC orbit but two free orbits: {pairs_same_lc_diff_free}"
          " (both > 0: neither partition refines the other)")
    # problem symmetry: fix L_0 and L_F. K_G = [[G,0],[G+G^-T, G^-T]].  Build via generators: CX(c,t) composed with a CZ/S correction.
    # Direct: all symplectic maps generated by the free gens that also fix L_F.  Enumerate the free group image (10,752) and filter.
    free_gens = [CX(a, b) for a in range(m) for b in range(m) if a != b] + [CZ(a, b) for a in range(m) for b in range(a + 1, m)] + [S(j) for j in range(m)]
    basis_imgs = tuple(1 << t for t in range(2 * m))
    group = {basis_imgs}; stack = [basis_imgs]
    while stack:
        M = stack.pop()
        for g in free_gens:
            N = tuple(g(v) for v in M)
            if N not in group:
                group.add(N); stack.append(N)
    stab = [M for M in group if frozenset(apply_sym(m, M, p) for p in lf) == lf]
    print(f"  free symplectic image {len(group)} (= 10,752); of these fixing L_F too: {len(stab)} (= |GL(3,2)| = 168)")
    orbs = {}
    for Z in Ls:
        key = min(tuple(sorted(apply_sym(m, M, p) for p in Z)) for M in stab)
        orbs.setdefault(key, []).append(Z)
    inv = {}
    for key, members in orbs.items():
        Z = members[0]
        r, e = beta_minus_dot_rank(m, Z)
        inv.setdefault((r, e), []).append(len(members))
    print(f"  orbits of the (L_0, L_F)-stabiliser on the 135 frames: {len(orbs)}; by (r, e): " +
          ", ".join(f"(r={r},e={e}): {sorted(v)}" for (r, e), v in sorted(inv.items())))
    # frames opposite L_0: Z = {(x, Gamma x)}; classify Gamma + I by (rank, alternating)
    top = [Z for Z in Ls if dist(m, l0, Z) == m]
    cls = {}
    mask = (1 << m) - 1
    for Z in top:
        G = [next(p & mask for p in Z if p >> m == (1 << j)) for j in range(m)]   # column j = Gamma e_j
        M = [G[j] ^ (1 << j) for j in range(m)]                                    # Gamma + I
        rk = rank_f2(M)
        alt = all((M[j] >> j & 1) == 0 for j in range(m))
        key = (rk, "alternating" if alt and rk else "non-alternating" if rk else "zero")
        cls.setdefault(key, []).append(Z)
    print(f"  frames opposite L_0 ((x|Gamma x), Gamma symmetric): {len(top)}; Gamma+I by (rank, type): " +
          ", ".join(f"{k}: {len(v)}" for k, v in sorted(cls.items())) +
          f"; d(Z, L_F) = rank(Gamma+I) for all: {all(dist(m, Z, lf) == k[0] for k, v in cls.items() for Z in v)}")
    # the stabiliser orbits on the top cell coincide with these classes?
    top_orbs = [members for members in orbs.values() if dist(m, l0, members[0]) == m]
    same = sorted(len(o) for o in top_orbs) == sorted(len(v) for v in cls.values())
    print(f"  stabiliser orbits on the top cell: sizes {sorted(len(o) for o in top_orbs)}; equal to the congruence classes of Gamma+I: {same}")

if __name__ == "__main__":
    check_A(int(sys.argv[1]) if len(sys.argv) > 1 else 9)
    check_B(4)
    check_C()
    print("DONE")
