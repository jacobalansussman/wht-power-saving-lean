# Order-2 sufficiency check (e3-slice6). Take an order-2 transfer triple and build the FULL slice program in which
# each move t_z is followed, on its own register, by the free element u_z = ((1-i) + (1+i) X^z)/2 of C[Sh]
# (one copy, one shift X^z, one addition), and the outputs are multiplied by v = sum_lambda i^|lambda| e_lambda
# (also in C[Sh]).  Evaluate at EVERY character psi of A: X^a -> (-1)^(lambda.a), lambda = psi mod 2.
import sys, itertools
import numpy as np
d = np.load(sys.argv[1])
A, B, Gam = d["A"], d["B"], d["Gam"]; dirs = [int(z) for z in d["dirs"]]; W = int(d["W"]); m = int(d["m"])
M = len(dirs)
dot = lambda z, y: bin(z & y).count("1") & 1
psis = set()
for ns in itertools.product(range(4), repeat=1 << m):
    psis.add(tuple(sum(n * dot(z, y) for n, y in zip(ns, range(1 << m))) % 4 for z in range(1 << m)))
worst = {"plain": 0.0, "corrected": 0.0}
for psi in sorted(psis):
    lam = sum((psi[1 << j] % 2) << j for j in range(m))
    target = (1j ** (sum(psi[1 << j] for j in range(m)) % 4)) * np.eye(W)
    w = np.array([1j ** psi[z] for z in dirs])
    sh = np.array([(-1) ** dot(lam, z) for z in dirs])         # X^z on the move's register
    u = ((1 - 1j) + (1 + 1j) * sh) / 2                           # the free C[Sh] element after each move
    v = 1j ** bin(lam).count("1")                                # sum_lambda i^|lambda| e_lambda at psi
    for name, c, scale in (("plain", w - 1, 1.0), ("corrected", w * u - 1, v)):
        out = scale * (np.eye(W) + A @ np.diag(c) @ np.linalg.solve(np.eye(M) - 0.5 * Gam @ np.diag(c), B))
        worst[name] = max(worst[name], np.abs(out - target).max())
print(f"characters {len(psis)}; max error without the C[Sh] corrections {worst['plain']:.3e}; "
      f"with them {worst['corrected']:.3e}")
