"""classify.py -- o3-shared (14 Oct 2026): minimal separating piece sets for W arrays with L = W + 1 end-pieces
over Q[R], R abelian (signs = translations Z2^k, order 3 x signs = Z3 x Z2^k, order 4 = Z4 (bit 4-cycle), Z4 x Z2).

Theory used (proved in NOTE.md section 1): Q[R] = product over Galois orbits O of characters of the fields Q(chi_O).
An L x W matrix P over Q[R] is, component by component, L vectors in Q(chi)^W.  Up to GL_W(Q[R]) (change of input
basis, Lemma 1) and units of Q[R] on the rows, a component with all L = W + 1 vectors spanning is determined by the
support S(O) of its unique linear relation (S = {k}: row k vanishes there; S = {i, j}: rows i, j parallel).
  separating            : P(chi) spans at every orbit (built in)
  minimal (no W rows form a basis)          <=> the S(O) have empty intersection
  rationed (no TWO rows extend to a basis)  <=> every pair {i, j} contains some S(O)       (Lemma 6 hypothesis;
                                               Theorem R allows at most one private end in 3-in-10 and 4-in-13)
Classes up to free equivalence = assignments O -> S(O) modulo row permutations S_L and Aut(R) (acting on orbits).
Exact check: explicit P = sum_O e_O P_O (central idempotents), regular representation ranks over Fractions.
usage: python3 -I classify.py"""
import itertools, sys, json
from fractions import Fraction
from math import gcd

# ------------------------------------------------------------------ abelian groups as Z_n1 x ... x Z_nk
def elements(ns):
    return list(itertools.product(*[range(n) for n in ns]))

def add(ns, a, b):
    return tuple((x + y) % n for x, y, n in zip(a, b, ns))

def automorphisms(ns):
    """all automorphisms as images of the standard generators (brute force, small groups)"""
    els = elements(ns)
    gens = [tuple(1 if i == j else 0 for i in range(len(ns))) for j in range(len(ns))]
    def order(a):
        k, x = 1, a
        while any(x):
            x = add(ns, x, a); k += 1
        return k
    cand = [[e for e in els if order(e) == ns[j] or (ns[j] % order(e) == 0 and False)] for j in range(len(ns))]
    cand = [[e for e in els if (order(e) == ns[j])] for j in range(len(ns))]
    auts = []
    for imgs in itertools.product(*cand):
        def phi(a):
            out = tuple(0 for _ in ns)
            for j, aj in enumerate(a):
                for _ in range(aj):
                    out = add(ns, out, imgs[j])
            return out
        im = {phi(a) for a in els}
        if len(im) == len(els):
            # check well defined: relations n_j * g_j = 0 hold since order(img) = n_j
            auts.append({a: phi(a) for a in els})
    return auts

def characters(ns):
    """characters as tuples of exponents k_j: chi(a) = exp(2 pi i sum_j k_j a_j / n_j)"""
    return elements(ns)

def char_order(ns, k):
    L = 1
    for kj, n in zip(k, ns):
        o = n // gcd(kj, n)
        L = L * o // gcd(L, o)
    return L

def galois_orbits(ns):
    N = 1
    for n in ns:
        N = N * n // gcd(N, n)
    units = [u for u in range(1, N + 1) if gcd(u, N) == 1]
    seen = set(); orbs = []
    for k in characters(ns):
        if k in seen:
            continue
        orb = sorted({tuple((u * kj) % n for kj, n in zip(k, ns)) for u in units})
        for x in orb:
            seen.add(x)
        orbs.append(orb)
    return orbs

def aut_on_orbits(ns, orbs):
    """permutations of orbit indices induced by Aut(R) (dual action on characters)"""
    els = elements(ns)
    idx = {}
    for i, o in enumerate(orbs):
        for k in o:
            idx[k] = i
    perms = set()
    for phi in automorphisms(ns):
        # chi o phi: exponent vector of the character a -> chi(phi(a)); compute by values on generators
        p = []
        for o in orbs:
            k = o[0]
            # chi(phi(e_j)) = exp(2 pi i sum_m k_m phi(e_j)_m / n_m); find k' with that value on generators
            N = 1
            for n in ns:
                N = N * n // gcd(N, n)
            def val(a):  # in Z_N
                return sum(km * am * (N // nm) for km, am, nm in zip(k, a, ns)) % N
            kp = []
            for j, n in enumerate(ns):
                ej = tuple(1 if i == j else 0 for i in range(len(ns)))
                v = val(phi[ej])          # = k'_j * N / n_j
                assert v % (N // n) == 0
                kp.append(v // (N // n))
            p.append(idx[tuple(kp)])
        perms.add(tuple(p))
    return sorted(perms)

# ------------------------------------------------------------------ combinatorial classification
def valid(assign, L, cond):
    if cond in ("min", "rat2"):
        inter = set(range(L))
        for S in assign:
            inter &= set(S)
        if inter:
            return False
    if cond == "rat2":
        for i, j in itertools.combinations(range(L), 2):
            if not any(set(S) <= {i, j} for S in assign):
                return False
    return True

def canon(assign, L, operms):
    best = None
    for rp in itertools.permutations(range(L)):
        ren = [tuple(sorted(rp[x] for x in S)) for S in assign]
        for op in operms:
            t = tuple(ren[op[i]] for i in range(len(assign)))   # orbit i receives the set of orbit op[i]
            if best is None or t < best:
                best = t
    return best

def classify(ns, W, cond="rat2", maxenum=3_000_000):
    L = W + 1
    orbs = galois_orbits(ns)
    operms = aut_on_orbits(ns, orbs)
    subsets = [S for r in range(1, L + 1) for S in itertools.combinations(range(L), r)]
    tot = len(subsets) ** len(orbs)
    if tot > maxenum:
        return {"orbits": len(orbs), "enumerated": False, "assignments": tot}
    classes = {}
    nvalid = 0
    for assign in itertools.product(subsets, repeat=len(orbs)):
        if not valid(assign, L, cond):
            continue
        nvalid += 1
        c = canon(assign, L, operms)
        classes[c] = classes.get(c, 0) + 1
    return {"orbits": len(orbs), "orbit_sizes": [len(o) for o in orbs], "enumerated": True,
            "valid": nvalid, "classes": len(classes), "reps": sorted(classes)}

def core_types(W):
    """minimal covers of all pairs of [L] by sets of size <= 2 (a singleton {k} covers every pair through k,
    a doubleton covers itself), up to S_L: the 'rationed core' each class must contain"""
    L = W + 1
    pairs = list(itertools.combinations(range(L), 2))
    small = [(k,) for k in range(L)] + pairs
    out = set()
    for r in range(1, L + 2):
        for fam in itertools.combinations(small, r):
            cov = all(any(set(S) <= set(p) for S in fam) for p in pairs)
            if not cov:
                continue
            if any(all(any(set(S) <= set(p) for S in fam[:i] + fam[i + 1:]) for p in pairs) for i in range(r)):
                continue      # not minimal
            out.add(min(tuple(sorted(tuple(sorted(rp[x] for x in S)) for S in fam)) for rp in itertools.permutations(range(L))))
    return sorted(out, key=lambda f: (len(f), f))

# ------------------------------------------------------------------ explicit matrices and exact checks
def rank_frac(rows):
    M = [list(r) for r in rows]; r = 0
    if not M: return 0
    nc = len(M[0])
    for c in range(nc):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        piv = M[r][c]
        M[r] = [x / piv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r

def explicit(ns, W, assign, orbs):
    """P over Q[R] with prescribed S(O) at each orbit, rational entries via the orbit idempotents
    e_O = (1/|R|) sum_a (sum_{chi in O} chi(a)^-1) a  (rational because O is a Galois orbit)."""
    import cmath
    L = W + 1
    els = elements(ns)
    def chi(k, a):
        return cmath.exp(2j * cmath.pi * sum(kj * aj / n for kj, aj, n in zip(k, a, ns)))
    idem = []
    for o in orbs:
        coeffs = []
        for a in els:
            v = sum(chi(k, a).conjugate() for k in o) / len(els)
            fr = Fraction(round(v.real * len(els) * 8), len(els) * 8)
            assert abs(v.imag) < 1e-9 and abs(float(fr) - v.real) < 1e-9
            coeffs.append(fr)
        idem.append(coeffs)
    def local(S):
        """L vectors in Q^W with relation support S (rational: same vectors at every character of the orbit)"""
        S = list(S); rest = [k for k in range(L) if k not in S]
        vec = [[0] * W for _ in range(L)]
        # circuit on S spanning coordinates 0..|S|-2, the rest a basis of the remaining coordinates
        for t, k in enumerate(S[:-1]):
            vec[k][t] = 1
        for t in range(len(S) - 1):
            vec[S[-1]][t] = 1 if len(S) > 1 else 0
        for t, k in enumerate(rest):
            vec[k][len(S) - 1 + t] = 1
        return vec
    P = [[[Fraction(0)] * len(els) for _ in range(W)] for _ in range(L)]
    for oi, S in enumerate(assign):
        vec = local(S)
        for k in range(L):
            for l in range(W):
                if vec[k][l]:
                    P[k][l] = [x + vec[k][l] * y for x, y in zip(P[k][l], idem[oi])]
    return P

def span_rows(ns, P_rows):
    els = elements(ns); ix = {a: i for i, a in enumerate(els)}
    out = []
    for row in P_rows:
        for g in els:
            v = []
            for ent in row:
                w = [Fraction(0)] * len(els)
                for i, c in enumerate(ent):
                    if c:
                        w[ix[add(ns, g, els[i])]] += c
                v += w
            out.append(v)
    return out

def check(ns, W, P):
    n = len(elements(ns)); L = len(P)
    full = rank_frac(span_rows(ns, P)) == W * n
    wsub = [rank_frac(span_rows(ns, [P[k] for k in S])) == W * n for S in itertools.combinations(range(L), W)]
    pair = [rank_frac(span_rows(ns, [P[i], P[j]])) == 2 * n for i, j in itertools.combinations(range(L), 2)]
    single = [rank_frac(span_rows(ns, [P[i]])) == n for i in range(L)]
    terms = [[sum(1 for c in ent if c) for ent in row] for row in P]
    dens = sorted({c.denominator for row in P for ent in row for c in ent if c})
    return {"separates": full, "some_W_subset_basis": any(wsub), "some_pair_extends": any(pair),
            "unimodular_rows": sum(single), "terms_per_entry": terms, "denominators": dens}

GROUPS = {
    "Z2 (one sign)": [2], "Z2^2 (signs)": [2, 2], "Z2^3 (signs)": [2, 2, 2],
    "Z3 (order 3 alone)": [3], "Z6 = Z3 x Z2 (order 3 x one sign)": [3, 2], "Z3 x Z2^2 (order 3 x signs)": [3, 2, 2],
    "Z4 (order 4)": [4], "Z4 x Z2 (order 4 x sign)": [4, 2],
}

if __name__ == "__main__":
    out = {}
    for W in (3, 4):
        print(f"=== W = {W}, L = {W + 1}", flush=True)
        ct = core_types(W)
        print(f"rationed cores (minimal covers of the pairs of rows by S-sets of size <= 2), up to row permutation: {len(ct)}")
        for c in ct:
            print("   ", [list(s) for s in c])
        out[f"cores_W{W}"] = [[list(s) for s in c] for c in ct]
        for name, ns in GROUPS.items():
            for cond in (("min", "rat2") if W == 3 else ("rat2",)):
                r = classify(ns, W, cond)
                key = f"W{W} {name} {cond}"
                if not r["enumerated"]:
                    print(f"  {name:38s} {cond:4s}: {r['orbits']} orbits, {r['assignments']} assignments (not enumerated)")
                    out[key] = r; continue
                print(f"  {name:38s} {cond:4s}: orbits {r['orbits']} sizes {r['orbit_sizes']}: valid assignments "
                      f"{r['valid']}, classes up to S_L x Aut(R): {r['classes']}", flush=True)
                out[key] = {k: v for k, v in r.items() if k != "reps"}
                if cond == "rat2" and 0 < r["classes"] <= 12:
                    orbs = galois_orbits(ns)
                    for rep in r["reps"]:
                        P = explicit(ns, W, rep, orbs)
                        ck = check(ns, W, P)
                        print(f"      class {[list(s) for s in rep]}: {ck}")
                        out.setdefault(f"W{W} {name} reps", []).append({"S": [list(s) for s in rep], "check": ck})
    json.dump(out, open("out_classify.json", "w"), indent=1, default=str)
