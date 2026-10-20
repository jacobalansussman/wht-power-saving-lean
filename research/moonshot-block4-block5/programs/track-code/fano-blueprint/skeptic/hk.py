"""v-ms-fano: re-check ms-fano's hook depth values with tb4-tight's blk9 logic (via ms-fano's blkk.py copy, k from
HOOKK) at a chosen PRIME, and for depth-3 elements in several different direction triples (uniformity check).
usage: PRIME=p HOOKK=k python3 -I hk.py dump NAME z1,z2,.. ; nullsp M_NAME.bin null_NAME.bin ; python3 -I hk.py load NAME"""
import sys, os, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import blkk as B
here = os.path.dirname(os.path.abspath(__file__))
mode, name = sys.argv[1], sys.argv[2]
if mode == 'dump':
    zs = tuple(int(z) for z in sys.argv[3].split(','))
    t0 = time.time(); M = B.block0(B.prod(zs))
    with open(os.path.join(here, f"M_{name}.bin"), 'wb') as fh:
        np.array([M.shape[0], M.shape[1], B.p], dtype=np.uint32).tofile(fh); M.astype(np.uint32).tofile(fh)
    print(f"p={B.p} k={B.k} {name} dirs {zs} block {M.shape} ({time.time()-t0:.0f}s)", flush=True)
else:
    raw = np.fromfile(os.path.join(here, f"null_{name}.bin"), dtype=np.uint32, count=2)
    nf, C = int(raw[0]), int(raw[1]); r = C - nf
    print(f"RESULT p={B.p} k={B.k} {name}: rank {r}, 16*rank = {16*r}", flush=True)
    os.remove(os.path.join(here, f"null_{name}.bin")); os.remove(os.path.join(here, f"M_{name}.bin"))
