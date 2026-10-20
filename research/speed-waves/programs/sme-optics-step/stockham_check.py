# as-1: check the Stockham stage formula used by colfft.h (columns as the vector dimension)
import numpy as np
def stockham(x, radices, sign=-1):
    N = x.shape[0]; Ns = 1
    for R in radices:
        y = np.empty_like(x); NR = N // R
        for j in range(NR):
            k = j % Ns
            M = np.array([[np.exp(sign*2j*np.pi*((i*(o*Ns+k)) % (Ns*R))/(Ns*R)) for i in range(R)] for o in range(R)])
            v = np.array([x[j + i*NR] for i in range(R)])
            out = M @ v
            base = (j//Ns)*Ns*R + k
            for o in range(R): y[base + o*Ns] = out[o]
        x = y; Ns *= R
    return x
rng = np.random.default_rng(1)
for N, rad in ((1024, [2,8,8,8]), (2048, [4,8,8,8]), (4096, [8,8,8,8]), (64, [8,8])):
    x = rng.standard_normal((N, 3)) + 1j*rng.standard_normal((N, 3))
    y = stockham(x, rad); r = np.fft.fft(x, axis=0)
    print(N, rad, np.abs(y-r).max()/np.abs(r).max())
