# verify56.py (a2-gct): independent check of the relaxed 5-in-14 scheme in End(V_56) found by relax2.py.
# Rebuild the paid-move outputs u_1..u_14 from the saved coefficients (unnormalised: plain products, so the scheme is
# literally r_i = sum beta e + sum gamma u_j, u_i = h r_i), then solve DIRECTLY for output coefficients in S_56
# (one d x d block per type, acting on the type rows) with K e_l = sum_l' c e_l' + sum_i c' u_i, by least squares on
# the full 56 x 280 matrices (no variable projection), and report the relative residual. Also: the same with M = 13
# outputs dropped one at a time (each paid move must matter), and the Jacobian rank of the output map is skipped.
import numpy as np, relax2, sys
n, h, k, tl = relax2.load([8], np.random.default_rng(12345))
W, M = 5, 14
g = relax2.Game(n, h, k, tl, W, M)
x = np.load(sys.argv[1] if len(sys.argv) > 1 else 'sol_8_W5_M14_target_s1.npy')
beta, gamma, Uh, nrm = g.forward(x)
# unnormalised version: u_i = h r_i with gamma rescaled (same module)
U = []
scale = []
for i in range(M):
    r = np.zeros((n, g.N))
    for l in range(W): r[:, l * n:(l + 1) * n] += beta[i][l]
    for j in range(i): r += gamma[i, j] @ (U[j] / scale[j])
    u = h @ r; U.append(u); scale.append(1.0)   # nrm folded: uhat_j = u_j/|u_j| exactly as in forward
    U[-1] = u; scale[-1] = np.linalg.norm(u)
gens = [g.E[l] for l in range(W)] + [U[i] / scale[i] for i in range(M)]
assert all(np.allclose(U[i] / scale[i], Uh[i]) for i in range(M))
def solve(gens, target):
    tot = 0.0; res = 0.0
    for d, slots in tl:
        s = slots[0]
        A = np.vstack([G[s, :] for G in gens])           # (len(gens) d) x N : the type rows of each generator
        B = target[s, :]                                  # d x N
        # B = sum_g C_g A_g with C_g arbitrary d x d  <=>  rows of B in rowspace(A)
        C, *_ = np.linalg.lstsq(A.T, B.T, rcond=None)
        R = B - (A.T @ C).T
        res += float((R * R).sum()); tot += float((B * B).sum())
    return np.sqrt(res / tot)
worst = max(solve(gens, k @ g.E[l]) for l in range(W))
print(f'5 in 14 in End(V_56): max relative residual over the 5 arrays = {worst:.3e}')
for drop in range(M):
    gg = [G for t, G in enumerate(gens) if t != W + drop]
    print(f'  without paid move {drop + 1:2d}: max residual {max(solve(gg, k @ g.E[l]) for l in range(W)):.3e}')
cond = [np.linalg.svd(np.vstack([G[s[0], :] for G in gens]), compute_uv=False) for d, s in tl]
print('smallest singular values per type of the generator rows:', [float('%.3e' % c[min(len(c), g.N) - 1]) for c in cond])
print('coefficient norm |x| =', np.linalg.norm(x), ' max |entry| =', np.abs(x).max())
