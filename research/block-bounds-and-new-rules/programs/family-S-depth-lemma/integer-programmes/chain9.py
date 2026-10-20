# b4-depth4: f = 9 (k = 8) hook at m = 4, block chi = 0 (all 16 blocks isomorphic, see hook4odd.py), second prime
# p = 1048573, blocked exact rank (fastrank.py).  Computes rank of K (dem/16), and the plain-chain intersection
# dims dim(N_0 cap g1^-1 N cap ... cap gt^-1 N) on the block, t = 1..4, g = d1, d1d2, d1d2d4, K.
import sys, os, time
sys.argv = [sys.argv[0], os.environ.get('KK', '8')]
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, 'hook4odd.py')).read().split("D = {z: dirmat(z)")[0]
src = src.replace("p = 33554393", "p = 1048573")
exec(src)
sys.path.insert(0, here)
from fastrank import rank_simple, _mm
assert I * I % p == p - 1
D = {z: dirmat(z) for z in (1, 2, 4, 8)}
g1 = D[1]; g2 = mm(g1, D[2]); g3 = mm(g2, D[4]); K = mm(g3, D[8])
t0 = time.time()
Ms = [block0(A) for A in (g1, g2, g3, K)]
print(f"blocks built ({time.time()-t0:.0f}s)", flush=True)
for name, M in zip(('d1', 'g2', 'g3', 'K'), Ms):
    t = time.time(); r = rank_simple(M, p)
    print(f"  rank of {name} on block chi=0 = {r}  x16 = {16*r}  ({time.time()-t:.0f}s)", flush=True)
dims = [len(Ncols)]
for t in range(1, 5):
    tt = time.time(); r = rank_simple(np.vstack(Ms[:t]), p); dims.append(len(Ncols) - r)
    print(f"  chain t={t}: dim of intersection on block = {dims[-1]}  ({time.time()-tt:.0f}s)", flush=True)
g = [16 * (dims[t - 1] - dims[t]) for t in range(1, 5)]
print(f"plain chain f=9: intersection dims x16 = {[16*d for d in dims]}; gains {g}; sum {sum(g)}", flush=True)
