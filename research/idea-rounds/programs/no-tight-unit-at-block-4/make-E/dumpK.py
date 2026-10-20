# tb4-tight (resumed 14 Oct): dump the chi = 0 block of K_r (out x N for E; N x out for EP) for nullsp.c, using blk9.py
# (dumpM.py imported blk9_v2, which no longer exists); then 'load' converts nullsp.c's output to .npy.
# usage: python3 -I dumpK.py dump r E|EP ; python3 -I dumpK.py load r E|EP
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
mode, r, which = sys.argv[1], int(sys.argv[2]), sys.argv[3]
import blk9 as B
p = B.p; tag = f"r{r}_p{p}"
os.chdir(os.path.dirname(os.path.abspath(__file__)))
if mode == 'dump':
    t0 = time.time()
    if which == 'E': M = B.block0(B.prod(B.KBASES[r]))
    else: M = B.block0(B.prod(B.KBASES[r]), rows=B.Ncols, cols=B.Outrows)
    with open(f"M{which}_{tag}.bin", 'wb') as fh:
        np.array([M.shape[0], M.shape[1], p], dtype=np.uint32).tofile(fh); M.astype(np.uint32).tofile(fh)
    print(f"dumped {which} {M.shape} ({time.time()-t0:.0f}s)", flush=True)
else:
    raw = np.fromfile(f"null{which}_{tag}.bin", dtype=np.uint32)
    nf, C = int(raw[0]), int(raw[1]); X = raw[2:].reshape(nf, C).astype(np.int64)
    with open(f"M{which}_{tag}.bin", 'rb') as fh:
        hdr = np.fromfile(fh, dtype=np.uint32, count=3); M = np.fromfile(fh, dtype=np.uint32).reshape(hdr[0], hdr[1]).astype(np.int64)
    bad = 0
    for s in range(0, nf, 64):
        bad += int(((M @ X[s:s+64].T) % p).any())
    print(f"{which}: nullity {nf}; check M x = 0 on ALL basis vectors: {bad == 0}", flush=True)
    out = ("E" if which == 'E' else "EP") + f"_{tag}.npy"
    np.save(out + ".tmp.npy", X); os.replace(out + ".tmp.npy", out)
    print("saved", out)
