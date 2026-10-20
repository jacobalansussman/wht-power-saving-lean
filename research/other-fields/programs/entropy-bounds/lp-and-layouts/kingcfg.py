# kingcfg.py -- the king-graph relaxation: psi(v) = (d(u1,v), d(u2,v), d(u3,v)) maps Gamma_m non-expansively into the
# king graph {0..m}^3 (l-infinity adjacency). A layout's image there depends only on x_jk = d(u_j,u_k) and
# z_kj = d(u_k, u_j') (u_j' the partner of u_j). Test random admissible configurations: cut value > 3m - 1?
import sys, random, itertools, time
import w3image
from w3screen import exact_value
m = int(sys.argv[1]); N = int(sys.argv[2]); seed = int(sys.argv[3]); rng = random.Random(seed)
kn, kE, kidx = w3image.king(m); t0 = time.time(); vals = []
def cfg_pairs(x, z):
    # x[(j,k)] symmetric, z[(k,j)] = d(u_k, partner of u_j)
    a = [tuple(0 if k == j else x[tuple(sorted((j, k)))] for k in range(3)) for j in range(3)]
    b = [tuple(m if k == j else z[(k, j)] for k in range(3)) for j in range(3)]
    return [(kidx[a[j]], kidx[b[j]]) for j in range(3)]
cnt = 0
while cnt < N:
    x = {(0, 1): rng.randint(1, m), (0, 2): rng.randint(1, m), (1, 2): rng.randint(1, m)}
    z = {}
    for j in range(3):
        for k in range(3):
            if k != j: z[(k, j)] = rng.randint(max(1, m - x[tuple(sorted((j, k)))]), m)
    if any(z[(j, k)] == m and z[(k, j)] == m for j, k in ((0, 1), (0, 2), (1, 2))): continue   # compatible: skip
    # triangle inequalities among the u's
    if not (x[(0,1)] <= x[(0,2)] + x[(1,2)] and x[(0,2)] <= x[(0,1)] + x[(1,2)] and x[(1,2)] <= x[(0,1)] + x[(0,2)]): continue
    cnt += 1
    P = cfg_pairs(x, z)
    v = exact_value(w3image.cut_value(kn, kE, P))
    vals.append(float(v))
    if v <= 3 * m - 1: print("KING FAIL", x, z, float(v), flush=True)
print(f"m={m}: {N} random non-compatible configurations; king-graph exact value min {min(vals):.4f}, max {max(vals):.4f}; "
      f"> 3m-1 on {sum(1 for v in vals if v > 3*m-1)} ({time.time()-t0:.0f}s)")
