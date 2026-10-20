"""t_dsize.py: how big the free corrections are. Lambda(gbar) = [g] D(g) for an engine letter g; D(g) is a free element
of C[N]. For each letter type, the number of group elements (i^j X^a Z^b) in D(g) = sum_k D_k(g) e_k, and per block
the number of Pauli terms in D_k(g). Letters: paid dir z (every weight), phase letters diag(i^((z.x + c) mod 2)),
shifts X^a, elementary address maps x_i += x_j, random address maps, and K."""
import sys, os, random, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from complete import *
def cnt(E, g):
    Ds = E.Dk(g); return len(E.D(g)), [int(np.sum(np.abs(D) > 1e-12)) for D in Ds]
for m in (2, 3, 4, 5):
    E = Eng(m); n = E.n; rng = random.Random(m)
    rows = []
    for w in range(1, m + 1): rows.append((f"dir z, |z|={w}", E.dirm((1 << w) - 1)))
    for w in range(0, m + 1):
        for c in (0, 1):
            z = (1 << w) - 1
            rows.append((f"phase i^((z.x+{c}) mod 2), |z|={w}", np.diag([1j ** ((pc(z & x) + c) & 1) for x in range(n)])))
    rows.append(("shift X^a, a = all ones", E.Xm[n - 1]))
    M = np.zeros((n, n), complex)
    for x in range(n): M[x ^ (((x >> 1) & 1)), x] = 1          # x_0 += x_1
    rows.append(("address map x_0 += x_1", M))
    while True:
        cols = [rng.randrange(1, n) for _ in range(m)]
        if E.S.rank_rows(cols) == m: break
    rows.append(("random address map", E.lift_P(E.S.addr(cols))))
    rows.append(("K", E.K))
    for name, g in rows:
        tot, per = cnt(E, g)
        print(f"m={m} {name:34s}: group elements in D = {tot:5d}; Pauli terms in D_0..D_3 = {per}", flush=True)
