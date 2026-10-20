import numpy as np
exec(open('sym2core.py').read())
np.set_printoptions(linewidth=220, precision=4, suppress=True)
print('E6 denominators', sorted(set(x.denominator for row in E6q for (r, i) in row for x in (r, i))))
print('Ch =\n', np.round(Ch, 4)); print("R' =\n", np.round(Rp, 4)); print("K' =\n", np.round(Kp, 4))
def rk(X): return np.linalg.matrix_rank(X, tol=1e-9)
U = np.eye(10)[:6]; Up = np.eye(10)[6:]
print('rank R_UU, R_UU\', R_U\'U, R_U\'U\':', rk(Rp[:6, :6]), rk(Rp[:6, 6:]), rk(Rp[6:, :6]), rk(Rp[6:, 6:]))
print('dim U+UK', rk(np.concatenate([U, U @ Kp])), ' dim U\'+U\'K', rk(np.concatenate([Up, Up @ Kp])))
print('dim U+UR', rk(np.concatenate([U, U @ Rp])), ' dim U\'+U\'R', rk(np.concatenate([Up, Up @ Rp])))
S2 = Sys2(2, 3); rng = np.random.default_rng(1)
J = np.concatenate([S2.jac(rng.normal(size=S2.P) + 1j * rng.normal(size=S2.P)) for _ in range(2)], 1)
print('Sym2 problem 2 in 3: unknowns', S2.P, 'raw eqs', len(S2.Tv), 'dim T', rk(J), 'generic rank', rk(S2.jac(rng.normal(size=S2.P) + 1j * rng.normal(size=S2.P))))
