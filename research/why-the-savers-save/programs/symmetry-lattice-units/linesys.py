"""Line systems of the Clifford / Barnes-Wall / code family, as integer vectors (any common scaling).
Each builder returns an int64 array of representatives, one per line (up to sign)."""
import itertools, numpy as np

def _canon(rows):
    out = {}
    for r in rows:
        r = tuple(r); i = next(k for k, t in enumerate(r) if t)
        if r[i] < 0: r = tuple(-t for t in r)
        out[r] = 1
    return np.array(sorted(out), dtype=np.int64)

def affine_flats(n, d):
    """all d-dimensional affine subspaces of F_2^n, as sorted tuples of points."""
    pts = range(1 << n); seen = set(); out = []
    for basis in itertools.combinations(range(1, 1 << n), d):
        span = {0}
        for b in basis: span |= {s ^ b for s in span}
        if len(span) != 1 << d: continue
        for t in pts:
            F = tuple(sorted(s ^ t for s in span))
            if F not in seen: seen.add(F); out.append(F)
    return out

def quad_signs(F, n):
    """the sign patterns (-1)^{quadratic function} on the affine flat F (up to global sign)."""
    d = len(F).bit_length() - 1
    # coordinates on F: F = t + span(basis); recover by gaussian enumeration order
    pts = list(F); t = pts[0]; vecs = [p ^ t for p in pts]
    basis = []
    for v_ in vecs:
        x = v_
        for b in basis: x = min(x, x ^ b)
        if x: basis.append(v_)
    coord = {}
    for c in range(1 << d):
        p = 0
        for i in range(d):
            if c >> i & 1: p ^= basis[i]
        coord[p ^ t] = c
    mons = [()] + [(i,) for i in range(d)] + list(itertools.combinations(range(d), 2))
    pats = set()
    for coeffs in itertools.product((0, 1), repeat=len(mons)):
        if coeffs[0]: continue                      # drop the constant (global sign)
        s = []
        for p in pts:
            c = coord[p]; val = 0
            for a, m in zip(coeffs, mons):
                if a and all(c >> i & 1 for i in m): val ^= 1
            s.append(-1 if val else 1)
        pats.add(tuple(s))
    return pats

def stabilizer_lines(n, codims):
    """real stabilizer states of n qubits (vectors in Z^(2^n), unnormalised +-1 on an affine flat) whose support has
    codimension in `codims`, scaled to equal norm (support 2^d gets factor 2^((D-d)/2), D = max support dim used)."""
    N = 1 << n; rows = []
    dims = [n - c for c in codims]; D = max(dims)
    for d in dims:
        assert (D - d) % 2 == 0, 'mixed parities need sqrt(2) scaling'
        f = 2 ** ((D - d) // 2)
        for F in affine_flats(n, d):
            for s in quad_signs(F, n):
                x = [0] * N
                for p, sg in zip(F, s): x[p] = f * sg
                rows.append(x)
    return _canon(rows)

def bw16(): return stabilizer_lines(4, (1, 3))        # supports 8 and 2: 1920 + 240 = 2160 lines (= BW16 minimal vectors)
def bw16_other(): return stabilizer_lines(4, (0, 2, 4))  # supports 16, 4, 1: 1024 + 1120 + 16 (isometric to BW16)
def e8(): return stabilizer_lines(3, (1, 3))           # supports 4 and 1: 112 + 8 = 120 lines (E8 roots)
def bw32(): return stabilizer_lines(5, (1, 3, 5))      # supports 16, 4, 1: 73440 lines (BW32 minimal vectors)
def rm14(): return stabilizer_lines(4, (1,))           # support 8 only: the 30 RM(1,4) octads with 64 sign classes
def rm24w4(): return stabilizer_lines(4, (2,))         # support 4 only: the 140 weight-4 words of RM(2,4) x 8 signs

def golay_code():
    """extended binary Golay code [24,12,8] from the generator of the quadratic-residue construction."""
    q = [1, 3, 4, 5, 9, 11, 12, 13, 15, 16, 22]   # QRs mod 23 (nonzero squares) -> length-23 QR code, then extend
    qr = set((i * i) % 23 for i in range(1, 23))
    g = [0] * 23
    # generator polynomial approach is fiddly; use the idempotent-free description: rows = cyclic shifts of N(x) + ...
    words = set()
    base = [1 if i in qr or i == 0 else 0 for i in range(23)]          # 0 + QR positions
    rows = []
    for s in range(23):
        r = [base[(i - s) % 23] for i in range(23)]
        rows.append(r + [sum(r) % 2])
    rows.append([1] * 24)
    # span
    basis = []
    for r in rows:
        x = int(''.join(map(str, r)), 2)
        for b in basis: x = min(x, x ^ b)
        if x: basis.append(x)
    return basis

def leech():
    """Leech lattice minimal vectors (norm 32 in the standard sqrt(8)-scaled coordinates), 98280 lines."""
    basis = golay_code(); assert len(basis) == 12, len(basis)
    code = [0]
    for b in basis: code += [c ^ b for c in code]
    octads = [c for c in code if bin(c).count('1') == 8]; assert len(octads) == 759
    rows = []
    for c in octads:
        S = [i for i in range(24) if c >> (23 - i) & 1]
        for signs in itertools.product((1, -1), repeat=7):
            if signs.count(-1) % 2: continue
            x = [0] * 24; x[S[0]] = 2
            for sg, i in zip(signs, S[1:]): x[i] = 2 * sg
            rows.append(x)
    # the odd vectors (-3, 1^23) with signs: x = (1^24) - 4 e_i, then flip signs on a codeword's support
    for c in code:
        sg = [(-1 if c >> (23 - i) & 1 else 1) for i in range(24)]
        for i in range(24):
            x = [s for s in sg]; x[i] = -3 * sg[i]
            rows.append(x)
    for i, j in itertools.combinations(range(24), 2):
        for s in (4, -4):
            x = [0] * 24; x[i] = 4; x[j] = s; rows.append(x)
    return _canon(rows)
