"""hookk.py -- ms-fano: chi = 0 block ranks of the block-4 hook V (x) Lambda^k V for one element (product of dirs),
for any k (blkk.py = tb4-tight's blk9.py with k from HOOKK).  f = 16 * rank.  usage: HOOKK=k python3 -I hookk.py dump NAME z1,z2,.. ;
nullsp ; HOOKK=k python3 -I hookk.py load NAME"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import blkk as B
here = os.path.dirname(os.path.abspath(__file__))
mode, name = sys.argv[1], sys.argv[2]
if mode == 'dump':
    zs = tuple(int(z) for z in sys.argv[3].split(','))
    t0 = time.time()
    M = B.block0(B.prod(zs))
    with open(os.path.join(here, f"M_{name}.bin"), 'wb') as fh:
        np.array([M.shape[0], M.shape[1], B.p], dtype=np.uint32).tofile(fh); M.astype(np.uint32).tofile(fh)
    print(f"k={B.k} {name} dirs {zs} block {M.shape} ({time.time()-t0:.0f}s)", flush=True)
else:
    raw = np.fromfile(os.path.join(here, f"null_{name}.bin"), dtype=np.uint32, count=2)
    nf, C = int(raw[0]), int(raw[1]); r = C - nf
    print(f"RESULT k={B.k} {name}: rank {r}, f = {16*r}", flush=True)
    os.remove(os.path.join(here, f"null_{name}.bin")); os.remove(os.path.join(here, f"M_{name}.bin"))
