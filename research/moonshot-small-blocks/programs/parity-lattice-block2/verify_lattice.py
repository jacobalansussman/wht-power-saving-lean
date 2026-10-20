"""verify_lattice.py -- independent re-check (different code path: Hermite bases + membership by back substitution)
of the lattice potential claims at block 2, f = 2:
  Lam = O*1 + pi*Z{monomials of degree <= 3} + pi^2 L   is kept by every XW monomial generator, and
  lambda_Lam(D_z) = 6 (z = 01, 10, 11), lambda_Lam(K) = 10.
"""
import sys, itertools; sys.path.insert(0, '.')
import padic as P

def hermite(vecs, N):
    """lower-triangular-by-pivot basis (list of N vectors) of the O-span of vecs; pivots chosen with minimal
    valuation per coordinate, processed coordinate by coordinate."""
    rows = [list(v) for v in vecs]; basis = []
    for c in range(N):
        cand = [(P.vpi(r[c]), k) for k, r in enumerate(rows) if r[c] != P.ZERO]
        v, k = min(cand)
        p = rows.pop(k); pu = P.divpi(p[c], v)
        rows = [r if r[c] == P.ZERO else [P.gsub(P.gmul(pu, r[t]), P.gmul(P.divpi(r[c], v), p[t])) for t in range(N)] for r in rows]
        basis.append((c, v, p))
    return basis

def detval(basis): return sum(v for c, v, p in basis)

def member(basis, w, N):
    """is w in the O-span of the triangular basis? (back substitution; pivots in increasing coordinate order)"""
    w = list(w)
    for c, v, p in basis:
        x = w[c]
        if x == P.ZERO: continue
        vx = P.vpi(x)
        if vx < v: return False
        # coefficient a = x / p[c] = (x / pi^v) / unit ; multiply w by unit to avoid division: w <- pu*w - (x/pi^v) p
        pu = P.divpi(p[c], v); q = P.divpi(x, v)
        w = [P.gsub(P.gmul(pu, w[t]), P.gmul(q, p[t])) for t in range(N)]
    return all(t == P.ZERO for t in w)

def lattice_034(f):
    n = 2 * f; N = 1 << n; vecs = []
    vecs.append([P.ONE] * N)
    for d in range(0, 4):
        for S in itertools.combinations(range(n), d):
            vecs.append([P.PI if all(P.bit(x, b) for b in S) else P.ZERO for x in range(N)])
    for x in range(N):
        e = [P.ZERO] * N; e[x] = P.pipow(2); vecs.append(e)
    return vecs

def lam_det(s, G, basis, N):
    """length((pi^-s G Lam + Lam)/Lam) = v(det Lam) - v(det(pi^-s G Lam + Lam)) = v(det Lam) + N s - v(det(G Lam + pi^s Lam))"""
    cols = [p for c, v, p in basis]
    imgs = []
    for p in cols:
        imgs.append([sum_row(G[i], p) for i in range(N)])
    sc = [[P.gmul(P.pipow(s), t) for t in p] for p in cols]
    return detval(basis) + N * s - detval(hermite(imgs + sc, N))

def sum_row(row, vec):
    acc = P.ZERO
    for a, b in zip(row, vec):
        if a != P.ZERO and b != P.ZERO: acc = P.gadd(acc, P.gmul(a, b))
    return acc

def apply(m, vec):
    out = [P.ZERO] * len(vec)
    for x in range(len(vec)):
        out[m.sig[x]] = P.gmul(P.IPOW[m.ph[x]], vec[x])
    return out

if __name__ == '__main__':
    f = 2; N = 16
    basis = hermite(lattice_034(f), N)
    print("index of Lam in L: 2^%d  (expected 15 + 16 - 15 = ... computed from pivots)" % (detval(basis)))
    gens = P.xw_generators(f)
    allinv = True


    for name, m in sorted(gens.items()):
        ok = all(member(basis, apply(m, p), N) for c, v, p in basis)
        allinv &= ok
        if not ok: print("  NOT kept by", name)
    print("kept by all %d XW monomial generators:" % len(gens), allinv)
    for z in (1, 2, 3):
        s, G = P.paid(f, z)
        print("lambda(D_%d) =" % z, lam_det(s, G, basis, N))
    s, G = P.kernel(f)
    print("lambda(K) =", lam_det(s, G, basis, N))
    # the standard lattice for comparison
    bL = hermite([[P.ONE if y == x else P.ZERO for y in range(N)] for x in range(N)], N)
    print("standard L: lambda(D_1) =", lam_det(*P.paid(f, 1), bL, N), " lambda(K) =", lam_det(*P.kernel(f), bL, N))
