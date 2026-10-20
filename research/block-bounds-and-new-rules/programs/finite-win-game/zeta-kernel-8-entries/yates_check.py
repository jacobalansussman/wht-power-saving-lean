# sanity: the algebra solver finds the Yates word (12 calls, pi = id, all c = 1) and rejects a corrupted one
import sys; sys.path.insert(0, '.')
import algebra
N = 8
calls = [(S, S | (1 << k)) for k in range(3) for S in range(N) if not S >> k & 1]
pi = list(range(N))
eqs, nz = algebra.equations(N, pi, calls)
print("Yates calls", len(calls), "equations", len(eqs))
for p in (31, 37):
    s = algebra.solve_fp(len(calls), eqs, nz, p, want_all=True)
    print(" p", p, "solutions", len(s), "e.g.", s[:2])
print(" exact rational check with c=1:", algebra.check_word_exact(N, pi, calls, [1] * len(calls)))
bad = calls[:-1]
eqs, nz = algebra.equations(N, pi, bad)
print("Yates minus last call: equations", None if eqs is None else len(eqs), "(None = structurally impossible)")
# leaves of the G=12 DFS for pi=id: how many are solvable mod 31 (must be > 0, Yates is among them)
hits = n = 0
for line in open('leaves_b3_g12_id.txt'):
    n += 1
    if n > 3000: break
    pi2, c2 = algebra.parse_leaf(line)
    e, z = algebra.equations(N, pi2, c2)
    if e is not None and algebra.solve_fp(len(c2), e, z, 31): hits += 1
print("G=12 pi=id leaves checked", min(n, 3000), "solvable mod 31:", hits)
eqs, nz = algebra.equations(N, pi, calls)
print("binomial test on Yates (must be None = consistent):", algebra.binomial_inconsistent(len(calls), eqs))
# toy inconsistent system: c0*c1 = c2 (equality), c0*c1 + c2 = 0 (zero entry)  ->  2 c2 = 0, impossible
print("binomial test on toy inconsistent system (must give a relation):",
      algebra.binomial_inconsistent(3, [[(1, 0b011), (-1, 0b100)], [(1, 0b011), (1, 0b100)]]))
