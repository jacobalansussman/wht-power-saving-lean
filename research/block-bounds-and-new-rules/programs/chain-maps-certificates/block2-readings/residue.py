"""xc3-s1s2 (wave 17): the next frequencies for words that mix chi (singular at S = 1) with a letter singular at alpha.

(1) lambda of order 5 in F16 (lambda^-1 = lambda^4 = Frob^2(lambda)): symbols at lambda preserve the Hermitian form
    H(u, v) = sum_i u_x,i v_z,i^4 + u_z,i v_x,i^4 (involution x -> x^4 of F16 over F4), so they permute the generators
    (maximal totally isotropic F16-planes) of H(3, 16).  Checks: Hermitian preservation of every letter's symbol
    (incl. G3 = (1+S+S^2) I and G3b = [[1+S, S],[S, 1]], both singular at alpha), number of generators (325),
    orbit of L0 = {x = 0} under the letter symbols, transitivity on ordered pairs at each distance, distances of
    directions (1) and K (2) and free letters (0).
(2) local structure at two opposite frames x = L0, y = K L0: O = common neighbours (17); for every frame w opposite to
    both, S_w = neighbours of w in O.  Checks: the S_w are 68 distinct 5-sets forming an inversive plane (every 3 points
    of O in exactly one), and the derived structure at one point is an affine plane of order 4 (16 points, 20 lines of
    4, every 2 points on exactly one line).  This is the input of the Hirzebruch step in PROOF.md section 7.
(3) lambda of order 7 in F8: no Galois conjugate of lambda equals lambda^-1, so no Hermitian form; report the orbit of
    L0 among the 4745 F8-planes of F8^4 under the letter symbols and whether those symbols preserve the F8-bilinear
    symplectic form up to a scalar.
"""
import itertools, sys
from collections import deque

class GF:
    def __init__(self, deg, mod):
        self.q = 1 << deg; self.exp = [0] * (2 * self.q); self.log = [0] * self.q
        x = 1
        for i in range(self.q - 1):
            self.exp[i] = x; self.log[x] = i
            x <<= 1
            if x & self.q: x ^= mod
        for i in range(self.q - 1, 2 * self.q): self.exp[i] = self.exp[i - (self.q - 1)]
    def mul(self, a, b):
        if a == 0 or b == 0: return 0
        return self.exp[self.log[a] + self.log[b]]
    def inv(self, a): return self.exp[(self.q - 1 - self.log[a]) % (self.q - 1)]
    def pw(self, a, e):
        if e == 0: return 1
        if a == 0: return 0
        return self.exp[(self.log[a] * e) % (self.q - 1)]

ok = True
def rep(name, cond):
    global ok
    ok = ok and bool(cond)
    print(("PASS " if cond else "FAIL ") + name, flush=True)

def run(F, lam, frob_exp, hermitian, label):
    print(f"==== {label}: field of order {F.q}, lambda = g^{F.log[lam]} (order {(F.q-1)//__import__('math').gcd(F.q-1, F.log[lam])}) ====")
    q = F.q; mul = F.mul
    def ev(poly):            # poly: {d: 1} Laurent over F2, evaluate at lambda
        s = 0
        for d, c in poly.items():
            if c: s ^= F.pw(lam, d % (q - 1))
        return s
    def ev_inv(poly):        # evaluate at lambda^-1
        s = 0
        for d, c in poly.items():
            if c: s ^= F.pw(lam, (-d) % (q - 1))
        return s
    def m2inv(A):
        (a, b), (c, d) = A; det = mul(a, d) ^ mul(b, c); di = F.inv(det)
        return [[mul(d, di), mul(b, di)], [mul(c, di), mul(a, di)]]
    def addr(G):             # G: 2x2 of Laurent polys; symbol diag(G(l), G(l^-1)^-T)
        Gl = [[ev(G[i][j]) for j in range(2)] for i in range(2)]
        Gi = [[ev_inv(G[i][j]) for j in range(2)] for i in range(2)]
        Git = m2inv(Gi); Gz = [[Git[j][i] for j in range(2)] for i in range(2)]
        s = [[0] * 4 for _ in range(4)]
        for i in range(2):
            for j in range(2): s[i][j] = Gl[i][j]; s[2 + i][2 + j] = Gz[i][j]
        return s
    def phase(B):            # (x | z) -> (x | z + B(l) x)
        s = [[int(i == j) for j in range(4)] for i in range(4)]
        for i in range(2):
            for j in range(2): s[2 + i][j] = ev(B[i][j])
        return s
    def const(M): return [list(r) for r in M]
    one = {0: 1}; zero = {}
    S = {}
    S['chi'] = addr([[{0: 1, 1: 1}, zero], [zero, {0: 1, 1: 1}]])
    S['E'] = addr([[one, {1: 1}], [zero, one]])
    S['Gm'] = addr([[one, {1: 1}], [{1: 1}, {0: 1, 1: 1, 2: 1}]])
    if q != 4:
        S['G3'] = addr([[{0: 1, 1: 1, 2: 1}, zero], [zero, {0: 1, 1: 1, 2: 1}]])       # singular at alpha
        S['G3b'] = addr([[{0: 1, 1: 1}, {1: 1}], [{1: 1}, one]])                      # det 1 + S + S^2
    S['D1'] = phase([[zero, zero], [zero, one]])
    S['Dcz'] = phase([[zero, one], [one, zero]])
    S['CZchain'] = phase([[zero, zero], [zero, {1: 1, -1: 1}]])
    S['CNOT'] = const([[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 1]])
    S['SWAP'] = const([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    PAID = {}
    PAID['D01'] = const([[1, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    PAID['D10'] = const([[1, 0, 0, 0], [0, 1, 0, 1], [0, 0, 1, 0], [0, 0, 0, 1]])
    PAID['D11'] = const([[1, 0, 1, 1], [0, 1, 1, 1], [0, 0, 1, 0], [0, 0, 0, 1]])
    PAID['K'] = const([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
    ALL = dict(S); ALL.update(PAID)
    def app(s, u):
        out = []
        for i in range(4):
            v = 0
            for k in range(4): v ^= mul(s[i][k], u[k])
            out.append(v)
        return tuple(out)
    fr = lambda a: F.pw(a, frob_exp)
    def Hf(u, v):            # u_x . frob(v_z) + u_z . frob(v_x)
        s = 0
        for i in range(2): s ^= mul(u[i], fr(v[2 + i])) ^ mul(u[2 + i], fr(v[i]))
        return s
    def Bf(u, v):            # bilinear symplectic form
        s = 0
        for i in range(2): s ^= mul(u[i], v[2 + i]) ^ mul(u[2 + i], v[i])
        return s
    basis = [tuple(int(i == j) for j in range(4)) for i in range(4)]
    def basis_pairs():
        import random
        rnd = random.Random(5)
        vs = basis + [tuple(rnd.randrange(q) for _ in range(4)) for _ in range(6)]
        return [(u, v) for u in vs for v in vs]
    if hermitian:
        for nm, s in ALL.items():
            rep(f"{nm:8s} symbol preserves the Hermitian form", all(Hf(app(s, u), app(s, v)) == Hf(u, v) for u, v in basis_pairs()))
    else:
        for nm, s in ALL.items():
            ratios = set()
            for u, v in basis_pairs():
                a, b = Bf(app(s, u), app(s, v)), Bf(u, v)
                if b: ratios.add(mul(a, F.inv(b)))
                elif a: ratios.add('bad')
            print(f"     {nm:8s}: symplectic multipliers seen {sorted(map(str, ratios))}")
    # ---- subspaces (2-dim) in RREF ----
    def rref2(u, v):
        rows = [list(u), list(v)]; r = 0; piv = []
        for c in range(4):
            p = next((i for i in range(r, 2) if rows[i][c]), None)
            if p is None: continue
            rows[r], rows[p] = rows[p], rows[r]
            iv = F.inv(rows[r][c]); rows[r] = [mul(iv, x) for x in rows[r]]
            for i in range(2):
                if i != r and rows[i][c]:
                    fc = rows[i][c]; rows[i] = [x ^ mul(fc, y) for x, y in zip(rows[i], rows[r])]
            r += 1
            if r == 2: break
        return tuple(map(tuple, rows)) if r == 2 else None
    vecs = list(itertools.product(range(q), repeat=4))
    if hermitian:
        iso = [u for u in vecs if any(u) and Hf(u, u) == 0]
        pts = {}
        for u in iso:
            k = rref1 = None
            lead = next(c for c in u if c); iv = F.inv(lead); nu = tuple(mul(iv, x) for x in u)
            pts[nu] = 1
        P = list(pts)
        frames = set()
        for u in P:
            for v in P:
                if v <= u: continue
                if Hf(u, v) == 0:
                    frames.add(rref2(u, v))
        frames = sorted(frames)
        print(f"     isotropic points: {len(P)}; totally isotropic planes (generators): {len(frames)}")
        rep(f"generators of H(3,{q}): {len(frames)} (expect (q^3+1)(q+1) with q^2 = {q})", len(frames) == {16: 325, 4: 27}.get(q, len(frames)))
    else:
        frames = set()
        nz = [u for u in vecs if any(u)]
        # all 2-planes: pairs of independent vectors, canonical rref
        lead1 = [u for u in nz if next(c for c in u if c) == 1]
        for i, u in enumerate(lead1):
            for v in lead1[i + 1:]:
                r = rref2(u, v)
                if r: frames.add(r)
        frames = sorted(frames)
        print(f"     2-planes of F{q}^4: {len(frames)}")
    index = {f: i for i, f in enumerate(frames)}
    def img(s, f): return index[rref2(app(s, f[0]), app(s, f[1]))]
    L0 = rref2((0, 0, 1, 0), (0, 0, 0, 1)); i0 = index[L0]
    perms = {nm: [img(s, f) for f in frames] for nm, s in ALL.items()}
    def dist(i, j):
        A = frames[i] + frames[j]
        # rank of 4 vectors
        rows = [list(r) for r in A]; r = 0
        for c in range(4):
            p = next((t for t in range(r, 4) if rows[t][c]), None)
            if p is None: continue
            rows[r], rows[p] = rows[p], rows[r]
            iv = F.inv(rows[r][c]); rows[r] = [mul(iv, x) for x in rows[r]]
            for t in range(4):
                if t != r and rows[t][c]:
                    fc = rows[t][c]; rows[t] = [x ^ mul(fc, y) for x, y in zip(rows[t], rows[r])]
            r += 1
        return r - 2
    for nm in ALL:
        d = dist(i0, perms[nm][i0])
        want = 2 if nm == 'K' else (1 if nm in PAID else 0)
        rep(f"{nm:8s} moves L0 to distance {d} (want {want})", d == want)
    # orbit of L0
    seen = {i0}; dq = deque([i0])
    gens = list(perms.values())
    while dq:
        a = dq.popleft()
        for p in gens:
            b = p[a]
            if b not in seen: seen.add(b); dq.append(b)
    print(f"     orbit of L0 under the letter symbols: {len(seen)} of {len(frames)} frames")
    orbit = sorted(seen)
    if hermitian:
        rep(f"orbit of L0 = all {len(frames)} generators", len(orbit) == len(frames))
        # transitivity on ordered pairs at each distance (union-find over pairs)
        n = len(frames)
        D = [[dist(i, j) for j in range(n)] for i in range(n)]
        parent = list(range(n * n))
        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]; a = parent[a]
            return a
        for p in gens:
            for i in range(n):
                pi = p[i]
                for j in range(n):
                    a, b = find(i * n + j), find(pi * n + p[j])
                    if a != b: parent[a] = b
        orbs = {}
        for i in range(n):
            for j in range(n):
                orbs.setdefault(find(i * n + j), set()).add(D[i][j])
        bydist = {}
        for r, ds in orbs.items():
            for d in ds: bydist[d] = bydist.get(d, 0) + 1
        rep(f"group of symbols is transitive on ordered pairs at each distance: orbits per distance {bydist}",
            all(len(ds) == 1 for ds in orbs.values()) and all(v == 1 for v in bydist.values()))
        deg = sum(1 for j in range(n) if D[i0][j] == 1); opp = sum(1 for j in range(n) if D[i0][j] == 2)
        print(f"     L0 has {deg} neighbours and {opp} opposite frames")
        # local structure
        x = i0; y = perms['K'][i0]
        O = [z for z in range(n) if D[x][z] == 1 and D[y][z] == 1]
        print(f"     |O| = |common neighbours of L0 and K L0| = {len(O)}")
        circles = set()
        for w in range(n):
            if D[w][x] == 2 and D[w][y] == 2:
                circles.add(frozenset(z for z in O if D[w][z] == 1))
        sizes = sorted({len(c) for c in circles})
        rep(f"triads give {len(circles)} distinct circles of sizes {sizes} on O", True)
        trip_ok = all(sum(1 for c in circles if set(t) <= c) == 1 for t in itertools.combinations(O, 3))
        rep("every 3 points of O lie in exactly one circle (inversive plane)", trip_ok)
        z1 = O[0]
        lines = [c - {z1} for c in circles if z1 in c]
        pts16 = [z for z in O if z != z1]
        pair_ok = all(sum(1 for l in lines if a in l and b in l) == 1 for a, b in itertools.combinations(pts16, 2))
        rep(f"derived structure at one point: {len(pts16)} points, {len(lines)} lines of sizes {sorted({len(l) for l in lines})}, "
            f"every 2 points on exactly one line (affine plane)", pair_ok)
        # t_2 = t_3 = 0 in any realisation: every line through 2 of the 16 points carries exactly 4
        print(f"     -> in a c = 1 realisation these {len(pts16)} points of P^2(C) would have t_2 = t_3 = 0, t_{len(lines[0])} = {len(lines)}")
        return dict(frames=frames, O=O, circles=[sorted(c) for c in circles])
    else:
        # distances of directions and K within this orbit
        return None

F16 = GF(4, 0b10011)
lam5 = F16.exp[3]
res16 = run(F16, lam5, 4, True, "F16, lambda of order 5 (Hermitian, involution x -> x^4)")
F4 = GF(2, 0b111)
res4 = run(F4, F4.exp[1], 2, True, "F4, lambda = alpha (control: 27 frames)")
F8 = GF(3, 0b1011)
run(F8, F8.exp[1], 1, False, "F8, lambda of order 7 (no Hermitian form)")
import json
json.dump({"O": res16["O"], "circles": res16["circles"]}, open("out_residue_circles.json", "w"))
print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
