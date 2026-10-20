"""veval.py: independent residual of a saved point. Data (targets T_l, directions z_i in time order) are taken from
rm2-solve's inst.canonical; the forms, parities and F = A_D (I + H_D)^-1 B_D^T (direct solve) are my own.
usage: veval.py SPEC ORDER DIRS(frames|all) X.npy"""
import sys, itertools
import numpy as np
sys.path.insert(0, '..')
from inst import canonical, slice_I
spec, order, dk, xf = sys.argv[1:5]
m = 5
kind, rest = spec.split(':')
if kind == 'k':
    ks = [int(v) for v in rest.split(',')]
    inst, fd = canonical(m, ks, dirs=dk, order=order, seed=0)
else:
    inst, fd = slice_I(m, int(rest), dirs=dk, order=order, seed=0)
T = [np.array(t) for t in inst.T]; dirs = inst.dirs; W, M = len(T), len(dirs)
# check targets = multiplication matrices by alpha^(2k) in their basis
F = fd.F
if kind == 'k':
    for k, Tl in zip(ks, T):
        Mk = np.array([[F.tr(F.mul(F.mul(fd.fb[i], F.alpha(2 * k)), fd.fb[j])) for j in range(m)] for i in range(m)])
        assert (Mk == Tl).all(), "target mismatch"
x = np.load(xf); nA = W * M
A = x[:nA].reshape(W, M); B = x[nA:2 * nA].reshape(W, M)
H = np.zeros((M, M), dtype=x.dtype); H[np.triu_indices(M, 1)] = x[2 * nA:]
Z = np.array([[(z >> j) & 1 for j in range(m)] for z in dirs])
tot = 0.0; per = np.zeros(W); nf = 0; worst = 0.0
# every symmetric matrix S over F2 is a form Q(z) = z^T S z restricted... use Q(z) = sum_{i<=j} c_ij z_i z_j
pairs = [(i, j) for i in range(m) for j in range(i, m)]
for c in itertools.product([0, 1], repeat=len(pairs)):
    if not any(c): continue
    t = np.array([sum(ci * z[i] * z[j] for ci, (i, j) in zip(c, pairs)) % 2 for z in Z])
    p = np.array([sum(ci * Tl[i, j] for ci, (i, j) in zip(c, pairs)) % 2 for Tl in T])
    Dm = np.nonzero(t)[0]
    if len(Dm) == 0:
        assert not p.any(); continue
    nf += 1
    Fq = A[:, Dm] @ np.linalg.solve(np.eye(len(Dm)) + H[np.ix_(Dm, Dm)], B[:, Dm].T)
    R = Fq - np.diag(p)
    a2 = np.abs(R) ** 2
    per += a2.sum(1); tot += a2.sum(); worst = max(worst, np.abs(R).max())
print(f"{spec} {order} {dk} {xf}: W={W} M={M} forms={nf} residual {tot:.4e} max|entry| {worst:.3e}")
print("per-array", " ".join(f"{v:.3e}" for v in per))
