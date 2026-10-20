"""geo4.py -- ms-fano: block-4 Fano gadget frames in engine letters: rank-1 interval terms N<v> (15), and over one
shared rank-1 frame N<x> the seven rank-2 interval frames N<x,v> (v in F_2^4/x).  Prints the counts and timings."""
import sys, time
from fgeo import *
t0 = time.time()
m = 4
geo = Geo(m, seed=3); G = geo.G; D1 = G.D1
b1 = {}
tries = 0
while len(b1) < 15 and tries < 20000:
    tries += 1
    b = geo.rand_free(); g = G.mul(D1, b)
    if geo.d0(g) == 1 and geo.dK(g) == 3:
        s = geo.s_of(g)
        b1.setdefault(s, b)
print(f"rank-1 interval directions: {len(b1)} found after {tries} letters ({time.time()-t0:.1f}s)")
x = (1,)
gx = G.mul(D1, b1[x])
b2 = {}
tries = 0
while len(b2) < 7 and tries < 40000:
    tries += 1
    b = geo.rand_free(); g = G.mul(D1, b, gx)
    if geo.d0(g) == 2 and geo.dK(g) == 2:
        s = geo.s_of(g)
        if 1 in s:
            b2.setdefault(s, b)
print(f"over N<1>: rank-2 interval frames N<1,v>: {sorted(b2)} after {tries} letters ({time.time()-t0:.1f}s)")
