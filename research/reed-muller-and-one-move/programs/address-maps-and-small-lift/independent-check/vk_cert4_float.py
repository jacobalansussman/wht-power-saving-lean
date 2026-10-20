"""float (complex SVD) check that row-r caps of wedge4's cert equal the mod-p values (mod-p rank <= true rank, so
equality in float means the file's exact caps are reproduced as true ranks). usage: vk_cert4_float.py ROW cert_file"""
import itertools, sys, re
import numpy as np
row = int(sys.argv[1]); txt = open(sys.argv[2]).read()
blocks = re.findall(r"planes: (\[\[.*?\]\])\s*\n\s*cap for directions 1..15: (\[[^\]]*\])", txt)
planes = eval(blocks[row-1][0]); filecaps = eval(blocks[row-1][1])
n = 16; quads = list(itertools.combinations(range(n), 4)); qidx = {s: k for k, s in enumerate(quads)}
def sgn(a):
    s = 1
    for i in range(len(a)):
        for j in range(i+1, len(a)):
            if a[i] > a[j]: s = -s
    return s
def bvec(P):
    v = np.zeros(len(quads), dtype=complex); seen = set()
    for t in range(n):
        c = frozenset(t ^ x for x in P)
        if c in seen: continue
        seen.add(c); seq = [t ^ x for x in P]; v[qidx[tuple(sorted(seq))]] += sgn(seq)
    return v
N = np.array([bvec(P) for P in planes]).T
def dirm(w):
    M = np.zeros((n, n), dtype=complex)
    for x in range(n): M[x, x] += (1+1j)/2; M[x ^ w, x] += (1-1j)/2
    return M
Q = np.array(quads)
PERMS = list(itertools.permutations(range(4))); SG = [sgn(p) for p in PERMS]
def lam4(g):  # full 1820 x 1820 is big; apply to N's support only
    supp = sorted(set(np.nonzero(np.abs(N).sum(1) > 0)[0]))
    S = np.array([quads[k] for k in supp])
    out = np.zeros((len(quads), len(supp)), dtype=complex)
    for s, pm in zip(SG, PERMS):
        term = np.ones((len(quads), len(supp)), dtype=complex)
        for i in range(4): term = term * g[Q[:, i][:, None], S[:, pm[i]][None, :]]
        out += s * term
    return out @ N[supp, :]
def rk(A):
    sv = np.linalg.svd(A, compute_uv=False); return int((sv > 1e-9 * sv[0]).sum()), sv
dN, _ = rk(N)
caps = []; gaps = []
for w in range(1, n):
    r, sv = rk(np.hstack([N, lam4(dirm(w))])); caps.append(r - dN)
    gaps.append(sv[r-1] / max(sv[r] if r < len(sv) else 1e-300, 1e-300))
print(f"row {row}: dim N {dN}; float caps {caps}; file {filecaps}; equal {caps == filecaps}; min singular-value gap ratio {min(gaps):.2e}")
K = np.eye(n, dtype=complex)
for w in (1, 2, 4, 8): K = K @ dirm(w)
r, sv = rk(np.hstack([N, lam4(K)])); print(f"standard demand float: {r - dN}")
