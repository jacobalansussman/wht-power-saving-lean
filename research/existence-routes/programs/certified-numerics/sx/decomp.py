# decomp.py -- structure of A = span rho(B) on V(x)V = Sym^2 V + Lambda^2 V; does A split as A|Sym (+) A|Alt, and how do
# R = rho(D0) and Kr = rho(D1 D0) sit relative to it.
import numpy as np
exec(open('sxcore.py').read())
SW = np.zeros((16, 16))
for a in range(4):
    for b in range(4): SW[b * 4 + a, a * 4 + b] = 1
Ps = (np.eye(16) + SW) / 2; Pl = (np.eye(16) - SW) / 2
us, ss, _ = np.linalg.svd(Ps); Us = us[:, :10]; ul, sl, _ = np.linalg.svd(Pl); Ul = ul[:, :6]
def rk(X): s = np.linalg.svd(X, compute_uv=False); return int((s > 1e-9 * s[0]).sum())
As = np.array([(Us.conj().T @ M @ Us).ravel() for M in Pa]); Al = np.array([(Ul.conj().T @ M @ Ul).ravel() for M in Pa])
print('dim A|Sym2 =', rk(As), ' dim A|Alt2 =', rk(Al), ' dim A =', rk(Pa.reshape(dimA, -1)))
print('R, Kr commute with swap:', np.allclose(SW @ R, R @ SW), np.allclose(SW @ Kr, Kr @ SW))
# centre of A|Sym2: matrices in A|Sym2 commuting with all of A|Sym2
Bs = As.reshape(-1, 10, 10)
rows = []
for M in Bs:
    rows.append(np.kron(M, np.eye(10)) - np.kron(np.eye(10), M.T))  # vec(MX - XM) row-major
C = np.concatenate(rows, 0); _, sv, vh = np.linalg.svd(C); null = vh[sv.size - (sv < 1e-9 * sv[0]).sum():] if False else vh[(sv > 1e-9 * sv[0]).sum():]
print('commutant of A|Sym2 has dim', null.shape[0])
Z = null[0].reshape(10, 10) if null.shape[0] else None
if null.shape[0] >= 2:
    X = sum(np.random.randn() * v.reshape(10, 10) for v in null)
    ev = np.linalg.eigvals(X); print('eigenvalues of a random commutant element (multiplicities = irreducible dims):', np.round(np.sort_complex(ev), 4))
# Burnside on Alt2
print('A|Alt2 is all of End(Lambda^2) (36):', rk(Al) == 36)
Kl = Ul.conj().T @ Kr @ Ul
print('Kr restricted to Lambda^2 lies in A|Alt2:', rk(np.concatenate([Al, Kl.ravel()[None]], 0)) == rk(Al))
