"""fbound.py -- ml-logic (meta3, 14 Oct 2026). How "correct for EVERY f" becomes a finite check, engine game.

Column-uniform letters at block m generate a matrix group G in GL(V), V = C^(2^m); a register content is an element
c of C[G] and acts at f columns as  sum_g c_g g^(x)f.  Claim (proved in NOTE.md, Theorem F1):
  c acts as 0 at every f >= 1   <=>   c acts as 0 at f = 1..r   <=>   c = 0 in C[G]   (when G is finite),
where r = number of distinct values of tr(g) on G (Burnside-Brauer).  So in the engine game the f-bound is uniform
in the block (does not depend on the word).  This script:
  A. builds the engine group at m = 1, 2 exactly (letters: D_z = a I + b X^z (paid move, a=(1+i)/2, b=(1-i)/2),
     Paulis X_j, Z_j, phase S_j, CZ_jk, CNOT_jk, scalar i)  and counts distinct traces r(m);
  B. at m = 1 computes, for F = 1, 2, ..., the rank of the evaluation map C[G] -> (+)_{f<=F} End(V^(x)f) through its
     Gram matrix  Gram_F[g,h] = sum_{f=1..F} tr(g^-1 h)^f  (Hilbert-Schmidt), finds the first F* where the rank is |G|,
     and checks F* <= r; it then exhibits a content that is right at f = 1..F*-1 and wrong at f = F* (sharpness);
  C. at m = 3 samples random products to collect trace values (lower bound on r(3)) against the bound 1 + 8(2m+1).
Usage: python3 -I fbound.py      (small: < 200 MB, < 1 min)
"""
import sys
import numpy as np, itertools, random

SC = 1 << 12
def key(M):
    return (np.round(M.real * SC).astype(np.int64).tobytes() + np.round(M.imag * SC).astype(np.int64).tobytes())

def letters(m):
    n = 1 << m; I = np.eye(n, dtype=complex)
    a, b = (1 + 1j) / 2, (1 - 1j) / 2
    def X(z):
        M = np.zeros((n, n), complex)
        for x in range(n): M[x ^ z, x] = 1
        return M
    def diag(f): return np.diag([f(x) for x in range(n)]).astype(complex)
    L = {}
    for j in range(m):
        z = 1 << j
        L[f"D{j}"] = a * I + b * X(z)
        L[f"X{j}"] = X(z)
        L[f"Z{j}"] = diag(lambda x: (-1) ** ((x >> j) & 1))
        L[f"S{j}"] = diag(lambda x: 1j ** ((x >> j) & 1))
    for j, k in itertools.permutations(range(m), 2):
        if j < k: L[f"CZ{j}{k}"] = diag(lambda x: (-1) ** (((x >> j) & 1) * ((x >> k) & 1)))
        M = np.zeros((n, n), complex)
        for x in range(n): M[x ^ ((((x >> j) & 1)) << k), x] = 1
        L[f"CX{j}{k}"] = M
    L["iI"] = 1j * I
    return L

def bfs_group(m, cap=200000):
    L = letters(m); gens = list(L.values())
    n = 1 << m; I = np.eye(n, dtype=complex)
    elems = [I]; seen = {key(I): 0}; frontier = [I]
    while frontier:
        nxt = []
        for g in frontier:
            for s in gens:
                h = s @ g
                k = key(h)
                if k not in seen:
                    seen[k] = len(elems); elems.append(h); nxt.append(h)
                    if len(elems) > cap: return None
        frontier = nxt
    return elems

def distinct(vals, tol=1e-9):
    out = []
    for v in vals:
        if all(abs(v - w) > tol for w in out): out.append(v)
    return out

def partA_B():
    res = {}
    for m in (1, 2):
        G = bfs_group(m)
        tr = [np.trace(g) for g in G]
        r = len(distinct(np.round(np.array(tr), 9)))
        mods = sorted({round(abs(t) ** 2, 6) for t in tr})
        print(f"m = {m}: |G| = {len(G)}, distinct traces r = {r} (bound 1 + 8(2m+1) = {1 + 8 * (2 * m + 1)}), "
              f"|tr|^2 values {mods}", flush=True)
        res[m] = (G, r)
    G, r = res[1]
    N = len(G)
    idx = {key(g): i for i, g in enumerate(G)}
    inv = [idx[key(np.linalg.inv(g))] for g in G]
    # chi(g^-1 h) table
    Ttab = np.array([[np.trace(G[inv[i]] @ G[j]) for j in range(N)] for i in range(N)])
    Gram = np.zeros((N, N), complex); ranks = []
    Fstar = None; prevGram = None
    for F in range(1, r + 3):
        prevGram = Gram.copy()
        Gram = Gram + Ttab ** F
        ev = np.linalg.eigvalsh((Gram + Gram.conj().T) / 2)
        rk = int((ev > 1e-7 * ev.max()).sum())
        ranks.append(rk)
        if Fstar is None and rk == N: Fstar = F
    print(f"m = 1: rank of evaluation on C[G] for f = 1..F, F = 1..{r + 2}: {ranks}; |G| = {N}; first full F* = {Fstar}; "
          f"r = {r}; F* <= r: {Fstar is not None and Fstar <= r}", flush=True)
    # sharpness: a content c with c acting as 0 at f = 1..F*-1 but not at F*
    if Fstar and Fstar > 1:
        G0 = np.zeros((N, N), complex)
        for F in range(1, Fstar): G0 = G0 + Ttab ** F
        w, V = np.linalg.eigh((G0 + G0.conj().T) / 2)
        c = V[:, 0]       # smallest eigenvalue -> in the kernel of the f < F* evaluation
        def act(c, f):
            return sum(c[i] * np.kron_power(G[i], f) for i in range(N))
        def kron_power(M, f):
            out = np.ones((1, 1), complex)
            for _ in range(f): out = np.kron(out, M)
            return out
        np.kron_power = kron_power
        norms = [float(np.linalg.norm(act(c, f))) for f in range(1, Fstar + 1)]
        print(f"   sharpness: content c (|c| = 1) with ||sum c_g g^(x)f|| for f = 1..{Fstar}: "
              + ", ".join(f"{x:.2e}" for x in norms), flush=True)
    return res

def partC(m=3, samples=200000, seed=1):
    L = letters(m); gens = list(L.values()); rng = random.Random(seed)
    n = 1 << m; vals = set()
    g = np.eye(n, dtype=complex)
    for t in range(samples):
        g = gens[rng.randrange(len(gens))] @ g
        if t % 50 == 49: g = np.eye(n, dtype=complex)
        v = np.trace(g); vals.add((round(v.real, 6), round(v.imag, 6)))
    mods = sorted({round(x * x + y * y, 4) for x, y in vals})
    print(f"m = 3 (sampled {samples} random-walk elements): distinct traces seen {len(vals)} "
          f"(bound 1 + 8(2m+1) = {1 + 8 * (2 * m + 1)}); |tr|^2 values {mods}", flush=True)

if __name__ == "__main__":
    partA_B()
    partC()
