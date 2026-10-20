"""v-iso-quantum: independent small checks of iso-quantum's relations (numpy, < 100 MB, seconds).

(A) dir e1 = H S H exactly, and S dir S = omega H (omega = e^{i pi/4}): one paid direction = one Hadamard
    modulo free (Hadamard-free) letters, up to a global phase omega that is NOT in the engine's Z4 phases.
(B) For random Cliffords U on m = 2, 3 qubits: every column and every row of U has exactly 2^d nonzeros,
    d = frame distance d(L0, U L0). So the monomial rank of U is 2^d (Koenig on a 2^d-regular support),
    and of U^(x)f it is 2^(f d).  (relation 6's [derived here] identity)
(C) Frame distance is additive under tensor products: d(U (x) V) = d(U) + d(V). Hence the quantum H-count
    of K^(x)W (tensor copies) is exactly W m: no amortisation over tensor copies for Clifford targets.
    (my remark on relation 5)
(D) Support lemma behind relation 8's by-product, in a column-mixing game at m = 2, f = 1..4:
    letters = transversal shifts, Z4 phases, CS on a column pair (transversal), the chain map row ^= row << 1,
    the parity phase i^(parity of row 1); paid = transversal direction. Random words on 1 live + 2 scratch
    registers with random rational mixing. For every register content, nonzeros per column of its operator
    <= (number of paths) * 2^(f * max paid depth along a path). K^(x)f has 2^(m f) per column.
"""
import itertools
import numpy as np

rng = np.random.default_rng(7)
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
S = np.diag([1, 1j]).astype(complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
ca, cb = (1 + 1j) / 2, (1 - 1j) / 2
DIR = ca * I2 + cb * X
om = np.exp(1j * np.pi / 4)

print("(A) dir == H S H:", np.allclose(DIR, H @ S @ H), "; S dir S == omega H:", np.allclose(S @ DIR @ S, om * H))


def kron(*ms):
    out = np.array([[1.0 + 0j]])
    for M in ms:
        out = np.kron(out, M)
    return out


def one(m, j, g):
    return kron(*[g if k == j else I2 for k in range(m)])


def cnot(m, c, t):
    d = 2 ** m
    M = np.zeros((d, d), dtype=complex)
    for x in range(d):
        b = [(x >> (m - 1 - k)) & 1 for k in range(m)]
        if b[c]:
            b[t] ^= 1
        M[sum(v << (m - 1 - k) for k, v in enumerate(b)), x] = 1
    return M


def paulis(m):
    out = []
    for xs in itertools.product([0, 1], repeat=m):
        for zs in itertools.product([0, 1], repeat=m):
            out.append((xs + zs, kron(*[(X if a else I2) @ (Z if b else I2) for a, b in zip(xs, zs)])))
    return out


def lag_of_unitary(U, P):
    # the frame U L0 U^-1 = Paulis Q with U^-1 Q U diagonal (Z-type)
    L = []
    for key, Q in P:
        M = U.conj().T @ Q @ U
        if np.allclose(M, np.diag(np.diag(M))):
            L.append(key)
    return L


def rank2(vecs):
    basis = []
    for v in vecs:
        x = int("".join(map(str, v)), 2)
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x)
    return len(basis)


def frame_dist(U, m, P):
    L = lag_of_unitary(U, P)
    L0 = [k for k, _ in P if not any(k[:m])]
    inter = [k for k in L if k in set(L0) and any(k)]
    return m - rank2(inter)


def rand_clifford(m, length):
    gens = [one(m, j, H) for j in range(m)] + [one(m, j, S) for j in range(m)]
    gens += [cnot(m, a, b) for a in range(m) for b in range(m) if a != b]
    U = np.eye(2 ** m, dtype=complex)
    for _ in range(length):
        U = gens[rng.integers(len(gens))] @ U
    return U


ok_B = 0; tot_B = 0; hist = {}
for m in (2, 3):
    P = paulis(m)
    for _ in range(150):
        U = rand_clifford(m, int(rng.integers(0, 25)))
        d = frame_dist(U, m, P)
        nz = np.abs(U) > 1e-9
        cols, rows = set(nz.sum(axis=0)), set(nz.sum(axis=1))
        tot_B += 1
        if cols == {2 ** d} and rows == {2 ** d}:
            ok_B += 1
        hist[(m, d)] = hist.get((m, d), 0) + 1
print(f"(B) nonzeros per row and column = 2^d for {ok_B}/{tot_B} random Cliffords; (m, d) counts {sorted(hist.items())}")

ok_C = 0
P1, P2, P3 = paulis(1), paulis(2), paulis(3)
for _ in range(60):
    U = rand_clifford(1, int(rng.integers(0, 10)))
    V = rand_clifford(2, int(rng.integers(0, 20)))
    if frame_dist(np.kron(U, V), 3, P3) == frame_dist(U, 1, P1) + frame_dist(V, 2, P2):
        ok_C += 1
K2 = np.kron(DIR, DIR)
print(f"(C) d(U (x) V) = d(U) + d(V): {ok_C}/60 ; d(L0, K L0) at m = 2: {frame_dist(K2, 2, P2)}, "
      f"K (x) K (W = 2 tensor copies, 4 qubits) distance {frame_dist(np.kron(K2, K2), 4, paulis(4))}")

# (D) column-mixing support lemma, m = 2
m = 2


def addr_bits(x, f):
    # address x of m*f bits: column c holds bits (row 0, row 1); bit index = c*m + r, most significant first
    return [(x >> (m * f - 1 - i)) & 1 for i in range(m * f)]


def from_bits(b, f):
    return sum(v << (m * f - 1 - i) for i, v in enumerate(b))


def perm_matrix(fun, f):
    n = 2 ** (m * f)
    M = np.zeros((n, n), dtype=complex)
    for x in range(n):
        M[from_bits(fun(addr_bits(x, f)), f), x] = 1
    return M


def diag_matrix(fun, f):
    n = 2 ** (m * f)
    return np.diag([fun(addr_bits(x, f)) for x in range(n)]).astype(complex)


def letters(f):
    free = []
    # transversal shifts X^z
    for z in [(1, 0), (0, 1), (1, 1)]:
        free.append(perm_matrix(lambda b, z=z: [b[i] ^ z[i % m] for i in range(m * f)], f))
    # transversal CNOT (row1 ^= row0)
    free.append(perm_matrix(lambda b: [b[i] ^ (b[i - 1] if i % m == 1 else 0) for i in range(m * f)], f))
    # transversal S on row 0, CS on (row0,row1) per column (i^(sum over columns))
    free.append(diag_matrix(lambda b: 1j ** sum(b[c * m] for c in range(f)), f))
    free.append(diag_matrix(lambda b: 1j ** sum(b[c * m] * b[c * m + 1] for c in range(f)), f))
    # chain map row ^= row << 1 : column c gets column c-1 added (stops at column 0)
    def chain(b):
        out = list(b)
        for c in range(1, f):
            for r in range(m):
                out[c * m + r] ^= b[(c - 1) * m + r]
        return out
    free.append(perm_matrix(chain, f))
    # parity phase i^(parity of row 1 over all columns)
    free.append(diag_matrix(lambda b: 1j ** (sum(b[c * m + 1] for c in range(f)) % 2), f))
    dirs = []
    for z in [(1, 0), (0, 1), (1, 1)]:
        sh = perm_matrix(lambda b, z=z: [b[i] ^ z[i % m] for i in range(m * f)], f)
        # transversal direction dir z in every column = prod over columns of (ca + cb shift_z on that column)
        D = np.eye(2 ** (m * f), dtype=complex)
        for c in range(f):
            shc = perm_matrix(lambda b, z=z, c=c: [b[i] ^ (z[i % m] if i // m == c else 0) for i in range(m * f)], f)
            D = (ca * np.eye(2 ** (m * f)) + cb * shc) @ D
        dirs.append(D)
    return free, dirs


viol = 0; checked = 0; maxratio = 0.0
for f in (1, 2, 3, 4):
    free, dirs = letters(f)
    n = 2 ** (m * f)
    Kf = np.eye(1, dtype=complex)
    for _ in range(f):
        Kf = np.kron(Kf, np.kron(DIR, DIR))
    knz = int((np.abs(Kf) > 1e-9).sum(axis=0).max())
    for trial in range(40):
        # registers: (operator from live input, number of paths, max paid depth)
        regs = [(np.eye(n, dtype=complex), 1, 0), (np.zeros((n, n), dtype=complex), 0, 0), (np.zeros((n, n), dtype=complex), 0, 0)]
        for step in range(int(rng.integers(3, 14))):
            kind = rng.integers(3)
            r = int(rng.integers(3))
            op, npth, dep = regs[r]
            if kind == 0:
                regs[r] = (free[rng.integers(len(free))] @ op, npth, dep)
            elif kind == 1 and npth > 0:
                regs[r] = (dirs[rng.integers(len(dirs))] @ op, npth, dep + 1)
            else:
                # mixing: register r := a*r + b*s (rational coefficients)
                s = int(rng.integers(3))
                a, b = rng.integers(-2, 3, size=2)
                ops, nps, des = regs[s]
                if s == r:
                    regs[r] = ((a + b) * op, npth, dep)
                else:
                    regs[r] = (a * op + b * ops, npth * (a != 0) + nps * (b != 0), max(dep if a else 0, des if b else 0))
            for (op2, np2, de2) in regs:
                if np2 == 0:
                    continue
                nzc = int((np.abs(op2) > 1e-9).sum(axis=0).max())
                bound = np2 * 2 ** (f * de2)
                checked += 1
                maxratio = max(maxratio, nzc / bound)
                if nzc > bound:
                    viol += 1
    print(f"(D) f={f}: K^(x)f nonzeros per column {knz} = 2^(m f) = {2 ** (m * f)}")
print(f"(D) support bound nnz/col <= paths * 2^(f * depth): {checked - viol}/{checked} register states pass; max ratio {maxratio:.3f}")
