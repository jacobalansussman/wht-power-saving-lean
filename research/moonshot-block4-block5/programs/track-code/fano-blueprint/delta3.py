"""delta3.py -- ms-fano: the f = 9 hook value of a depth-3 element at block 4 (delta_3), to decide whether two
depth-3 targets fit in 5 rank-1 moves (2 delta_3 <= 5*13 in units of 2112?).  Uses tb4-tight's blk9.block0 (chi = 0
block, out reps x N reps) and nullsp.c; calibration: D1 -> rank 1716 (13*2112/16), D1D2 -> 3036 (23*2112/16).
usage: python3 -I delta3.py dump NAME z1,z2,..   (writes M_<NAME>.bin here) ; then nullsp ; then  load NAME"""
import sys, os, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../common'))
import numpy as np
import blk9 as B
here = os.path.dirname(os.path.abspath(__file__))
mode, name = sys.argv[1], sys.argv[2]
if mode == 'dump':
    zs = tuple(int(z) for z in sys.argv[3].split(','))
    t0 = time.time()
    M = B.block0(B.prod(zs))
    with open(os.path.join(here, f"M_{name}.bin"), 'wb') as fh:
        np.array([M.shape[0], M.shape[1], B.p], dtype=np.uint32).tofile(fh); M.astype(np.uint32).tofile(fh)
    print(f"dumped {name} dirs {zs} shape {M.shape} ({time.time()-t0:.0f}s)", flush=True)
else:
    raw = np.fromfile(os.path.join(here, f"null_{name}.bin"), dtype=np.uint32, count=2)
    nf, C = int(raw[0]), int(raw[1])
    r = C - nf
    print(f"{name}: chi=0 block {C} cols, nullity {nf}, rank {r}, f = 16*rank = {16*r} = {16*r/2112:.4f} x 2112", flush=True)
    os.remove(os.path.join(here, f"null_{name}.bin")); os.remove(os.path.join(here, f"M_{name}.bin"))
