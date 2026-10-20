"""v-xc3-s1s2-compute: label matrices from the actual unitaries (f = 4, 5; 8 and 10 qubits), and Lemma 4
(|tr U|^2 in {0, 2^dim ker(M(U) - I)}) on random words, incl. P = i^parity, chi, E, Gm, boundary H/S, k, l."""
import sys, random
src = open("vcheck.py").read().split('print("== (A) regularity')[0]
exec(src)
random.seed(11)
def qubit_of(j, r): return 2 * j + r
def unitary(name, f):
    n = 2 * f; N = 1 << n
    if name in ("chi", "E", "Gm", "k", "l", "CNOT", "SWAPQ", "chi_ac"):
        A = {"chi": lambda: addr_poly({(0, 0): [0, 1], (1, 1): [0, 1]}, f),
             "chi_ac": lambda: addr_poly({(0, 0): [0, -1], (1, 1): [0, -1]}, f),
             "E": lambda: addr_poly({(0, 0): [0], (0, 1): [1], (1, 1): [0]}, f),
             "Gm": lambda: addr_poly({(0, 0): [0], (0, 1): [1], (1, 0): [1], (1, 1): [0, 1, 2]}, f),
             "k": lambda: A_k(f), "l": lambda: A_l(f),
             "CNOT": lambda: addr_poly({(0, 0): [0], (1, 0): [0], (1, 1): [0]}, f),
             "SWAPQ": lambda: addr_poly({(0, 1): [0], (1, 0): [0]}, f)}[name]()
        U = np.zeros((N, N), complex)
        for s in range(N):
            b = np.array([(s >> p) & 1 for p in range(n)], np.int64)
            t = A.astype(np.int64) @ b % 2
            U[sum(int(t[p]) << p for p in range(n)), s] = 1
        return U
    if name in ("D1", "Dcz", "CZch", "P", "Slast"):
        pairs = {"D1": [(qubit_of(j, r), qubit_of(j, r)) for j in range(f) for r in range(2)],
                 "Dcz": [(qubit_of(j, 0), qubit_of(j, 1)) for j in range(f)],
                 "CZch": [(qubit_of(j, 0), qubit_of(j + 1, 1)) for j in range(f - 1)],
                 "P": [(qubit_of(j, 1), qubit_of(k, 1)) for j in range(f) for k in range(j, f)],
                 "Slast": [(qubit_of(f - 1, 0), qubit_of(f - 1, 0))]}[name]
        d = np.ones(N, complex)
        for s in range(N):
            e = 0
            for a, b in pairs:
                ba = (s >> a) & 1; bb = (s >> b) & 1
                e += (1 if a == b else 2) * ba * bb
            d[s] = 1j ** (e % 4)
        return np.diag(d)
    Hm = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    def on_qubits(qs):
        U = np.array([[1.0 + 0j]])
        for p in reversed(range(n)):
            U = np.kron(U, Hm if p in qs else np.eye(2))
        return U
    if name == "Hlast": return on_qubits({qubit_of(f - 1, 0)})
    if name == "Hfirst": return on_qubits({qubit_of(0, 1)})
    if name == "H0": return on_qubits({qubit_of(j, 0) for j in range(f)})
    if name == "H1": return on_qubits({qubit_of(j, 1) for j in range(f)})
    raise KeyError(name)
def pauli(n, xv, zv):
    X = np.array([[0, 1], [1, 0]], complex); Z = np.diag([1, -1]).astype(complex)
    U = np.array([[1.0 + 0j]])
    for p in reversed(range(n)):
        g = np.eye(2, dtype=complex)
        if xv[p]: g = g @ X
        if zv[p]: g = Z @ g if not xv[p] else Z @ X
        U = np.kron(U, g)
    return U
def label_of_unitary(U, f):
    n = 2 * f
    M = np.zeros((4 * f, 4 * f), np.uint8)
    for j in range(f):
        for t in range(4):
            xv = [0] * n; zv = [0] * n
            r = t % 2; p = qubit_of(j, r)
            (xv if t < 2 else zv)[p] = 1
            C = U @ pauli(n, xv, zv) @ U.conj().T
            found = False
            for col in range(4 * f):
                pass
            # find the Pauli: brute force over labels is too big; read x-part from support, z-part from phases
            nz = np.nonzero(np.abs(C[:, 0]) > 1e-9)[0]
            assert len(nz) == 1
            xs = int(nz[0])
            xbits = [(xs >> q) & 1 for q in range(n)]
            P0 = pauli(n, xbits, [0] * n)
            D = P0.conj().T @ C            # should be +- i^? Z^zbits (diagonal)
            dg = np.diag(D)
            assert np.allclose(np.abs(dg), 1) and np.allclose(D, np.diag(dg), atol=1e-9)
            zbits = [0] * n
            for q in range(n):
                zbits[q] = int(np.real(dg[1 << q] / dg[0]) < 0)
            for q in range(n):
                jj, rr = q // 2, q % 2
                M[lidx_x(jj, rr), 4 * j + t] = xbits[q]
                M[lidx_z(jj, rr), 4 * j + t] = zbits[q]
    return M
ok = True
for f in (4, 5):
    names = ["chi", "chi_ac", "E", "Gm", "k", "l", "CNOT", "SWAPQ", "D1", "Dcz", "CZch", "P", "Slast", "Hlast", "Hfirst", "H0", "H1"]
    US = {nm: unitary(nm, f) for nm in names}
    for nm in names:
        Mu = label_of_unitary(US[nm], f); Ml = LETTERS[nm][1](f)
        good = np.array_equal(Mu, Ml)
        ok &= good
        if not good: print(f"  f={f} {nm}: label of the unitary != label matrix used in vcheck.py")
    print(f"f={f}: label matrices of all {len(names)} letters agree with the unitaries' conjugation action: {ok}")
    bad = 0; stats = []
    allnames = names + [nm + "^-1" for nm in names]
    for _ in range(40):
        w = [random.choice(allnames) for _ in range(random.randint(1, 7))]
        U = np.eye(1 << (2 * f), dtype=complex); M = np.eye(4 * f, dtype=np.uint8)
        for nm in w:
            b = nm.replace("^-1", "")
            Ug = US[b] if not nm.endswith("^-1") else US[b].conj().T
            Mg = LETTERS[b][1](f) if not nm.endswith("^-1") else f2inv(LETTERS[b][1](f))
            U = U @ Ug; M = f2mul(M, Mg)
        t2 = abs(np.trace(U)) ** 2
        dk = 4 * f - f2rank(M ^ np.eye(4 * f, dtype=np.uint8))
        good = abs(t2) < 1e-6 or abs(t2 - 2 ** dk) < 1e-6 * 2 ** dk
        if not good: bad += 1; print("   LEMMA 4 FAIL", w, t2, 2 ** dk)
        stats.append((round(t2, 3), 2 ** dk))
    ok &= bad == 0
    print(f"f={f}: Lemma 4 on 40 random words: failures {bad}; sample (|tr|^2, 2^dimker):", stats[:10])
print("ALL PASS" if ok else "SOME FAIL")
