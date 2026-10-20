import numpy as np, warnings
warnings.filterwarnings("ignore")
from inst import canonical
from glm import Model, random_start
rng = np.random.default_rng(7)
inst, fd = canonical(3, [0, 1, 3], dirs="all", order="rand3")
md = Model(inst.t, inst.p, inst.W, inst.M); md.w = np.array([0.5, 2.0, 1.3])
x = random_start(md, rng, 1.0, 0.5)
def resid(x):
    F = md.Fval(x); ar = np.arange(md.W); F[:, ar, ar] -= inst.p; return F
r0 = resid(x); h = 1e-7
wr = np.broadcast_to(md.w[None, :, None], r0.shape).ravel()
J = np.zeros((r0.size, md.nvar), dtype=complex)
for v in range(md.nvar):
    e = np.zeros(md.nvar, dtype=complex); e[v] = h
    J[:, v] = ((resid(x + e) - resid(x - e)) / (2 * h)).ravel()
N, g, c = md.normal(x, chunk=11)
print("cost", c, (wr * np.abs(r0.ravel()) ** 2).sum(), md.cost(x))
print("N err", np.abs(N - J.conj().T @ (wr[:, None] * J)).max(), np.abs(N).max())
print("g err", np.abs(g - J.conj().T @ (wr * r0.ravel())).max(), np.abs(g).max())
