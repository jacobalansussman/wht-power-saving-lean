import time
from b5eng import *
t0=time.time()
P, ok = plain_chain_test(5, 2)
print("closure ok", ok, "dims", [(f[0], f[1]) for f in P.fin], "%.1fs" % (time.time()-t0))
al = np.ones(P.nal, dtype=complex)
for d in range(5, 0, -1):
    r, _ = P.vp(al, d, [], want_J=False)
    print("d=%d res^2=%.3e" % (d, np.sum(np.abs(r)**2)))
