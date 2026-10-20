"""t_check.py: check glm.py's F, gradient and J^H J against direct formulas (block 3, small sizes)."""
import numpy as np
from inst import canonical, slice_I
from glm import Model, random_start
rng = np.random.default_rng(5)
inst, fd = canonical(3, [0, 1], dirs="all", order="rand3")
md = Model(inst.t, inst.p, inst.W, inst.M)
x = random_start(md, rng, 1.0, 0.5)
A, B, H = md.unpack(x)
# direct F
err = 0
Fd = []
for q in range(inst.Nf):
    D = np.nonzero(inst.t[q])[0]
    Fq = A[:, D] @ np.linalg.inv(np.eye(len(D)) + H[np.ix_(D, D)]) @ B[:, D].T
    Fd.append(Fq)
Fd = np.array(Fd)
Fm = md.Fval(x)
print("F max err", np.abs(Fd - Fm).max())
# explicit Jacobian by complex finite differences (holomorphic): J[:, v] = (F(x + h e_v) - F(x))/h
def resid(x):
    F = md.Fval(x); ar = np.arange(md.W); F[:, ar, ar] -= inst.p; return F.ravel()
r0 = resid(x); h = 1e-7
J = np.zeros((r0.size, md.nvar), dtype=complex)
for v in range(md.nvar):
    e = np.zeros(md.nvar, dtype=complex); e[v] = h
    J[:, v] = (resid(x + e) - resid(x - e)) / (2 * h)
N, g, c = md.normal(x, chunk=7)
print("cost", c, (np.abs(r0) ** 2).sum())
print("N err", np.abs(N - J.conj().T @ J).max(), "scale", np.abs(N).max())
print("g err", np.abs(g - J.conj().T @ r0).max(), "scale", np.abs(g).max())
# plain solution 1 in 3 at block 3: A = B = 1, H = 2
inst1, fd = canonical(3, [0], dirs="frames", order="lex")
m1 = Model(inst1.t, inst1.p, 1, 3)
A = np.ones((1, 3), complex); B = np.ones((1, 3), complex); H = np.triu(2 * np.ones((3, 3)), 1).astype(complex)
print("plain 1 in 3 cost", m1.cost(m1.pack(A, B, H)), inst1.exact_tests())
