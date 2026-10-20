"""vfp.py -- v-ex2-modp: exhaustive F_p point count (alphas enumerated, betas by linear algebra mod p) for small systems."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlib import *
M = os.environ.get('MODP_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
def solvable(A, r, p):
    R = [row[:] + [b] for row, b in zip(A, r)]; n = len(A[0]) if A else 0; row = 0
    for c in range(n):
        piv = next((i for i in range(row, len(R)) if R[i][c] % p), None)
        if piv is None: continue
        R[row], R[piv] = R[piv], R[row]; inv = pow(R[row][c], -1, p)
        R[row] = [x * inv % p for x in R[row]]
        for i in range(len(R)):
            if i != row and R[i][c] % p:
                f = R[i][c]; R[i] = [(x - f * y) % p for x, y in zip(R[i], R[row])]
        row += 1
    return all(R[i][n] % p == 0 for i in range(row, len(R)))
def count(pre, p, maxn=10**6):
    names, F, _ = to_Q(read_ms(pre + '.p1.ms'), read_ms(pre + '.p2.ms'))
    ai = [i for i, v in enumerate(names) if v[0] == 'a']; bi = [i for i in range(len(names)) if i not in ai]
    Fp = [{k: (v.numerator * pow(v.denominator, -1, p)) % p for k, v in P.items()} for P in F]
    if p ** len(ai) > maxn: return None
    hits = 0
    for al in itertools.product(range(p), repeat=len(ai)):
        A = [[0] * len(bi) for _ in Fp]; r = [0] * len(Fp)
        for j, P in enumerate(Fp):
            for k, c in P.items():
                t = c
                for ii, a in zip(ai, al):
                    if k[ii]: t = t * pow(a, k[ii], p) % p
                bb = [q for q, b in enumerate(bi) if k[b]]
                assert len(bb) <= 1 and all(k[bi[q]] == 1 for q in bb)
                if bb: A[j][bb[0]] = (A[j][bb[0]] + t) % p
                else: r[j] = (r[j] - t) % p
        if solvable(A, r, p): hits += 1
    return hits, len(ai), len(bi)
for pre, ps in [(f'{M}/data/ms/t27_mid3_m2_s1.u', [2, 3]), (f'{M}/data/ms/cal_b3_2in5_s0.u', [2, 3, 5]),
                (f'{M}/data/ms/cal_b3_2in5_s2.u', [2, 3, 5]), (f'{M}/data/dec/t27_mid3_m2_s1/full_c_1_0_1_1.a0', [2, 3, 5, 7, 11, 13]),
                (f'{M}/data/dec/t27_mid3_m2_s1/full_c_1_0_1_0.a0', [2, 3, 5, 7])]:
    for p in ps:
        print(os.path.basename(pre), 'p =', p, '-> alphas with an F_p point (hits, n_alpha, n_beta):', count(pre, p), flush=True)
