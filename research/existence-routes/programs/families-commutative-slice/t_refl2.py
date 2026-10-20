import numpy as np, time, warnings
warnings.filterwarnings("ignore")
from refl import System
from refl2 import System2
rng = np.random.default_rng(3)
for (m, S, W) in [(3, [1,2,3,4,5,6,7], 2), (4, list(range(1, 16)), 4)]:
    s1 = System(m, S, W); s2 = System2(m, S, W)
    x = (rng.standard_normal(s1.nvar) + 1j*rng.standard_normal(s1.nvar))*0.5
    t = time.time(); r1, J1 = s1.residual(x, jac=True); t1 = time.time()-t
    t = time.time(); r2, J2 = s2.residual(x, jac=True); t2 = time.time()-t
    print(m, W, "res diff", np.abs(r1-r2).max(), "jac diff", np.abs(J1-J2).max(), "times", round(t1,3), round(t2,3), J2.shape)
