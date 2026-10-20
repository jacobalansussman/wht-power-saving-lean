# every binary quadratic form Q on F_2^m is a point evaluation of the f-column slice program for a suitable f:
# with columns y_1..y_f, the move t_z acts by i^{n(z)}, n(z) = #{c : z.y_c = 1}; choose two copies of e_j for each
# linear term and the triple e_j, e_k, e_j+e_k (plus two copies each of e_j, e_k) for each product term z_j z_k.
# Check: n(z) = 2 Q(z) mod 4 for every z (so the move acts by (-1)^Q(z)), shifts X^a = t_a^2 act trivially, and the
# kernel K = t_{e_1}...t_{e_m} acts by i^{sum_j n(e_j)} = (-1)^{tr Q}.
import itertools, random
def dot(a, b): return bin(a & b).count("1") & 1
for m in (3, 4, 5):
    bad = 0
    for trial in range(200):
        lin = [random.randint(0, 1) for _ in range(m)]
        quad = {(j, k): random.randint(0, 1) for j in range(m) for k in range(j + 1, m)}
        cols = []
        for j in range(m):
            if lin[j]: cols += [1 << j, 1 << j]
        for (j, k), q in quad.items():
            if q: cols += [1 << j, 1 << k, (1 << j) | (1 << k), 1 << j, 1 << j, 1 << k, 1 << k]
        Q = lambda z: (sum(lin[j] * ((z >> j) & 1) for j in range(m)) + sum(q * ((z >> j) & 1) * ((z >> k) & 1) for (j, k), q in quad.items())) % 2
        n = lambda z: sum(dot(z, y) for y in cols)
        ok = all(n(z) % 4 == 2 * Q(z) % 4 for z in range(1 << m))
        ok = ok and all((2 * n(a)) % 4 == 0 for a in range(1 << m))
        trQ = sum(Q(1 << j) for j in range(m)) % 2
        ok = ok and (sum(n(1 << j) for j in range(m)) % 4 == 2 * trQ % 4)
        bad += not ok
    print(f"m={m}: 200 random binary quadratic forms realised as point evaluations; failures {bad}")
