# Coefficient systems of call structures (leaves of tdfs).  For a leaf (pi, calls i<-j), every call is
# x_i <- x_i + c_t x_j with unknown c_t != 0 (WLOG, see the header of tdfs.c), and the word is correct iff
#   Pi(c) = P D Z  (exactly; no right diagonal is needed: a right diagonal is absorbed by conjugation),
# i.e. row i of Pi(c) vanishes off supp Z[pi(i)] and is CONSTANT (nonzero) on it.
# Entries of Pi(c) are sums of distinct squarefree monomials (one per path), all with coefficient +1.
# Exact decision over F_p: propagation (every equation is linear in each variable) + branching.
import sys, random, itertools

def target_rows(N):
    return [sum(1 << u for u in range(N) if (s & ~u) == 0) for s in range(N)]

def path_polys(N, calls):
    reg = [{i: [0]} for i in range(N)]          # reg[i][u] = list of monomials (bitmask over calls)
    for t, (i, j) in enumerate(calls):
        bit = 1 << t
        for u, ms in reg[j].items():
            reg[i].setdefault(u, []).extend(m | bit for m in ms)
    return reg

def equations(N, pi, calls):
    """returns (eqs, nonzero) ; each eq is a list of (coef, mono)"""
    ZR = target_rows(N)
    reg = path_polys(N, calls)
    eqs, nonzero = [], []
    for i in range(N):
        T = ZR[pi[i]]
        sup = [u for u in range(N) if T >> u & 1]
        for u, ms in reg[i].items():
            if not (T >> u & 1):
                eqs.append([(1, m) for m in ms])
        for u in sup:
            if u not in reg[i]:
                return None, None                  # structurally impossible
        u0 = sup[0]
        nonzero.append([(1, m) for m in reg[i][u0]])
        for u in sup[1:]:
            e = {}
            for m in reg[i][u]: e[m] = e.get(m, 0) + 1
            for m in reg[i][u0]: e[m] = e.get(m, 0) - 1
            e = [(c, m) for m, c in e.items() if c]
            if e: eqs.append(e)
    return eqs, nonzero

def evalpoly(poly, val, p):
    s = 0
    for c, m in poly:
        v = c
        k = 0
        while m:
            if m & 1: v = v * val[k] % p
            m >>= 1; k += 1
        s += v
    return s % p

def solve_fp(G, eqs, nonzero, p, want_all=False, limit=10**6):
    """exact search over F_p^*: returns list of solutions (each a tuple), [] if none."""
    sols = []
    allvars = (1 << G) - 1
    evars = []
    for e in eqs:
        vm = 0
        for c, m in e: vm |= m
        evars.append(vm)
    count = {v: sum(1 for vm in evars if vm >> v & 1) for v in range(G)}
    inv = [0] + [pow(a, p - 2, p) for a in range(1, p)]
    steps = [0]

    def propagate(val, assigned):
        changed = True
        while changed:
            changed = False
            for e, vm in zip(eqs, evars):
                un = vm & ~assigned
                if un == 0:
                    if evalpoly(e, val, p) != 0: return False
                    continue
                if un & (un - 1): continue
                v = un.bit_length() - 1
                a = b = 0
                for c, m in e:
                    t = c
                    mm = m & ~(1 << v); k = 0
                    while mm:
                        if mm & 1: t = t * val[k] % p
                        mm >>= 1; k += 1
                    if m >> v & 1: a += t
                    else: b += t
                a %= p; b %= p
                if a == 0:
                    if b != 0: return False
                    continue
                x = (-b) * inv[a] % p
                if x == 0: return False
                val[v] = x; assigned |= 1 << v; changed = True
        return assigned

    def rec(val, assigned):
        steps[0] += 1
        if steps[0] > limit: raise RuntimeError("limit")
        r = propagate(val, assigned)
        if r is False: return
        assigned = r
        if assigned == allvars:
            for nz in nonzero:
                if evalpoly(nz, val, p) == 0: return
            for e in eqs:
                if evalpoly(e, val, p) != 0: return
            sols.append(tuple(val))
            return
        v = max((v for v in range(G) if not assigned >> v & 1), key=lambda v: count[v])
        for x in range(1, p):
            nv = list(val); nv[v] = x
            rec(nv, assigned | (1 << v))
            if sols and not want_all: return
    rec([0] * G, 0)
    return sols

def parse_leaf(line):
    a, b = line.split("|")
    pi = [int(x) for x in a.split()[1:]]
    calls = [tuple(int(y) for y in x.split("<")) for x in b.split()[1:]]
    return pi, calls

def check_word_exact(N, pi, calls, cvals):
    """exact rational check: Pi(c) == P D Z, using Fractions"""
    from fractions import Fraction
    M = [[Fraction(int(i == j)) for j in range(N)] for i in range(N)]
    for (i, j), c in zip(calls, cvals):
        M[i] = [a + c * b for a, b in zip(M[i], M[j])]
    ZR = target_rows(N)
    for i in range(N):
        T = ZR[pi[i]]
        vals = set()
        for u in range(N):
            if T >> u & 1: vals.add(M[i][u])
            elif M[i][u] != 0: return False
        if len(vals) != 1 or 0 in vals: return False
    return True

# ---------------- exact binomial-lattice test (proves infeasibility over C when it fails) ----------------
def binomial_inconsistent(G, eqs):
    """Every 2-term equation reads m1 = s*m2 (s = +-1), i.e. c^alpha = s with alpha = e(m1)-e(m2).
    Over (C^*)^G these are jointly solvable iff every integer relation sum n_k alpha_k = 0 has prod s_k^n_k = 1.
    Returns a witness relation (list of (eq index, n_k)) if INconsistent, else None.  [exact integer algebra]"""
    rows, signs, idx = [], [], []
    for q, e in enumerate(eqs):
        if len(e) != 2: continue
        (c1, m1), (c2, m2) = e
        s = -c2 * c1            # c1 m1 + c2 m2 = 0, c's are +-1  ->  m1/m2 = -c2/c1 = -c2*c1
        alpha = [((m1 >> k) & 1) - ((m2 >> k) & 1) for k in range(G)]
        rows.append(alpha); signs.append(0 if s == 1 else 1); idx.append(q)
    k = len(rows)
    if k == 0: return None
    # integer row reduction of [A | I]; zero rows of A give the relation lattice
    M = [rows[i] + [1 if j == i else 0 for j in range(k)] for i in range(k)]
    r = 0
    for col in range(G):
        while True:
            piv = [i for i in range(r, k) if M[i][col] != 0]
            if not piv: break
            p = min(piv, key=lambda i: abs(M[i][col]))
            M[r], M[p] = M[p], M[r]
            done = True
            for i in range(r + 1, k):
                if M[i][col]:
                    f = M[i][col] // M[r][col]
                    M[i] = [a - f * b for a, b in zip(M[i], M[r])]
                    if M[i][col]: done = False
            if done: break
        if any(M[i][col] for i in range(r, k)):
            r += 1
    for i in range(r, k):
        if any(M[i][:G]): continue
        n = M[i][G:]
        par = sum(nk * sk for nk, sk in zip(n, signs)) % 2
        if par:
            return [(idx[j], n[j]) for j in range(k) if n[j]]
    return None

if __name__ == "__main__":
    N = int(sys.argv[1]); fn = sys.argv[2]
    primes = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [31, 37, 41, 43, 61, 73, 97, 101, 109, 113]
    n = struct_bad = binom_dead = hits = 0
    rest = []
    for line in open(fn):
        if not line.startswith("pi"): continue
        n += 1
        pi, calls = parse_leaf(line)
        eqs, nz = equations(N, pi, calls)
        if eqs is None:
            struct_bad += 1; continue
        w = binomial_inconsistent(len(calls), eqs)
        if w is not None:
            binom_dead += 1; continue
        found = []
        for p in primes:
            s = solve_fp(len(calls), eqs, nz, p)
            if s: found.append((p, s[0]))
        if found:
            hits += 1
            print("HIT", line.strip(), "primes", [p for p, _ in found], "first", found[0], flush=True)
        else:
            rest.append(line.strip())
    print(f"leaves {n}; structurally impossible {struct_bad}; killed EXACTLY by the binomial lattice test {binom_dead}; "
          f"with an F_p solution for some p in {primes}: {hits}; no F_p solution for any of these p (not an exact proof): {len(rest)}")
    for l in rest[:50]: print("  NOFP", l)
