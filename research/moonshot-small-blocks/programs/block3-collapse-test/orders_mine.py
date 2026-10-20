"""v-ms-xw3: independent orders of A.B (Clifford two-sided chain) and A.T2.B (with Toffoli), words applied left to
right (first letter first).  A: x_{r,j} ^= x_{r,j-1} (j >= 1, original values, every row r); B: x_{r,j} ^= x_{r,j+1};
T2: x_2 ^= x_0 x_1 in every column.  A.B is linear and acts on each row the same way, so its order is that of the
f x f GF(2) matrix; A.T2.B is computed as a permutation of the 8^f addresses."""
import sys, math
import numpy as np

def cycles_order(p):
    n = p.size; seen = np.zeros(n, bool); o = 1; ncyc = 0
    for s in range(n):
        if seen[s]: continue
        L = 0; x = s
        while not seen[x]:
            seen[x] = True; x = p[x]; L += 1
        o = o * L // math.gcd(o, L); ncyc += 1
    return o, ncyc

def perm_AB_row(f):          # one bit-row, 2^f states, bit j = digit column j
    x = np.arange(2 ** f, dtype=np.int64)
    bits = [(x >> j) & 1 for j in range(f)]
    a = [bits[0]] + [bits[j] ^ bits[j - 1] for j in range(1, f)]          # A
    b = [a[j] ^ a[j + 1] if j + 1 < f else a[j] for j in range(f)]        # then B
    return sum(b[j] << j for j in range(f))

def perm_AT2B(f):            # 8^f addresses, column j = bits 3j..3j+2, row r = bit r of the column
    x = np.arange(8 ** f, dtype=np.int64)
    c = [(x >> (3 * j)) & 7 for j in range(f)]
    a = [c[0]] + [c[j] ^ c[j - 1] for j in range(1, f)]                   # A on every row at once
    t = [d ^ (((d & 1) & ((d >> 1) & 1)) << 2) for d in a]                # T2
    b = [t[j] ^ t[j + 1] if j + 1 < f else t[j] for j in range(f)]        # B
    return sum(b[j] << (3 * j) for j in range(f))

fmaxAB = int(sys.argv[1]) if len(sys.argv) > 1 else 18
fmaxT = int(sys.argv[2]) if len(sys.argv) > 2 else 7
print("A.B (Clifford, two-sided chain): f, order, 2^f-1", flush=True)
for f in range(2, fmaxAB + 1):
    o, _ = cycles_order(perm_AB_row(f)); print(f"  f={f:2d} order {o:8d}   2^f-1 = {2**f-1}", flush=True)
print("A.T2.B (with Toffoli): f, order, cycles", flush=True)
for f in range(2, fmaxT + 1):
    o, nc = cycles_order(perm_AT2B(f)); print(f"  f={f:2d} order {o:8d} cycles {nc}", flush=True)
