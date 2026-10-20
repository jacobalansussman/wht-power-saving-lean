"""BW16 minimal vectors (4320, i.e. 2160 lines), Gram matrix (norm 4), pair classes, and F_2 ranks of the
orthogonality-pattern matrices.  numpy only.  Coordinates are indexed by F_2^4."""
import itertools, numpy as np

def affine_hyperplanes(n=4):
    pts = range(1 << n); out = []
    for a in range(1, 1 << n):
        for c in (0, 1):
            out.append(tuple(p for p in pts if (bin(a & p).count('1') & 1) == c))
    return out

def bw16_lines():
    L = []
    for i, j in itertools.combinations(range(16), 2):
        for s in (2, -2):
            x = [0]*16; x[i] = 2; x[j] = s; L.append(x)
    for H in affine_hyperplanes(4):
        for signs in itertools.product((1, -1), repeat=7):
            if signs.count(-1) % 2: continue
            x = [0]*16; x[H[0]] = 1
            for s, i in zip(signs, H[1:]): x[i] = s
            L.append(x)
    return np.array(L, dtype=np.int64)

def gram(L):
    return (L @ L.T) // 2

def f2_rank(M):
    """rank over F_2 of a 0/1 numpy matrix (bit-packed python ints)."""
    rows = [int(''.join('1' if b else '0' for b in r), 2) for r in (M & 1).astype(bool)]
    basis = {}
    for x in rows:
        while x:
            h = x.bit_length() - 1
            if h in basis: x ^= basis[h]
            else: basis[h] = x; break
    return len(basis)

if __name__ == '__main__':
    L = bw16_lines(); v = len(L)
    G = gram(L)
    print('lines', v, 'norms', set(np.diag(G)), 'rank_Q', np.linalg.matrix_rank(L.astype(float)))
    vals, cnt = np.unique(G[np.triu_indices(v, 1)], return_counts=True)
    print('off-diagonal ips (norm-4 scale):', dict(zip(vals.tolist(), cnt.tolist())))
    row0 = G[0]; print('ip profile of line 0:', dict(zip(*np.unique(row0, return_counts=True))))
    row1 = G[300]; print('ip profile of line 300:', dict(zip(*np.unique(row1, return_counts=True))))
    np.save('bw16_L.npy', L)
