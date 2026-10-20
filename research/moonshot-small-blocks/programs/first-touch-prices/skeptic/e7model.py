# Skeptic's independent model of Sp(6,2) = W(E7)^+ (the rotation subgroup of the Weyl group of E7).
# Nothing from c2-finite-group is imported.  The 7-dim representation: transvection t_v -> -s_alpha, alpha the root
# with alpha mod (2L + radical) = v.  Larger modules: tensor products of the 7.  Everything is integer.
import itertools, numpy as np
# --- E7 root system in simple-root coordinates (Bourbaki: 1-3-4-5-6-7 chain, 2 attached to 4) ---
EDGES = [(1, 3), (3, 4), (4, 5), (5, 6), (6, 7), (2, 4)]
C = 2 * np.eye(7, dtype=np.int64)
for a, b in EDGES: C[a - 1, b - 1] = C[b - 1, a - 1] = -1
def ip(x, y): return int(np.asarray(x) @ C @ np.asarray(y))
simple = [tuple(int(i == j) for j in range(7)) for i in range(7)]
roots = set(simple); fr = list(simple)
while fr:
    nx = []
    for r in fr:
        for s in simple:
            t = tuple(np.array(r) - ip(r, s) * np.array(s))
            if t not in roots: roots.add(t); nx.append(t)
    fr = nx
roots = sorted(roots); pos = [r for r in roots if all(c >= 0 for c in r)]
assert len(roots) == 126 and len(pos) == 63
def refl(a):
    a = np.array(a, dtype=np.int64); return np.eye(7, dtype=np.int64) - np.outer(a, C @ a)
# --- mod 2: F_2^7 with b(x, y) = x^T C y mod 2; radical; symplectic basis; labels (x|z) as 6-bit ints ---
def b2(x, y): return ip(x, y) & 1
vecs = [tuple((i >> j) & 1 for j in range(7)) for i in range(128)]
rad = [v for v in vecs if all(b2(v, w) == 0 for w in vecs)]
assert len(rad) == 2, len(rad)
def add(x, y): return tuple((a + b) & 1 for a, b in zip(x, y))
# symplectic Gram-Schmidt
E, F = [], []; pool = [v for v in vecs if any(v)]
def orth_ok(v): return all(b2(v, e) == 0 for e in E) and all(b2(v, f) == 0 for f in F)
for _ in range(3):
    cand = [v for v in pool if orth_ok(v) and v not in rad]
    e = cand[0]; f = next(w for w in cand if b2(e, w) == 1)
    E.append(e); F.append(f)
def label(v):   # (x|z): x_i = b(v, f_i), z_i = b(v, e_i); bits 0..2 = x, 3..5 = z
    v = tuple(int(c) & 1 for c in v)
    return sum(b2(v, F[i]) << i for i in range(3)) | sum(b2(v, E[i]) << (3 + i) for i in range(3))
def Bl(u, v):   # x.z' + x'.z on labels
    return (bin((u & 7) & (v >> 3)).count('1') + bin((v & 7) & (u >> 3)).count('1')) & 1
for v in vecs:
    for w in vecs: assert Bl(label(v), label(w)) == b2(v, w)
root_of = {}
for r in pos:
    l = label(r); assert l != 0 and l not in root_of; root_of[l] = r
assert len(root_of) == 63
def T(w):        # transvection on labels as a function
    return lambda v: v ^ (w if Bl(v, w) else 0)
# check: the reduction of s_alpha is t_label(alpha)
for l, r in root_of.items():
    S = refl(r); t = T(l)
    for v in vecs:
        assert label(S @ np.array(v)) == t(label(v))
def R7_t(w): return -refl(root_of[w])            # rho_7 of the transvection t_w (det +1)
def word_R7(ws):                                 # rho_7(t_w1 t_w2 ... t_wr), acting on columns
    M = np.eye(7, dtype=np.int64)
    for w in ws: M = M @ R7_t(w)
    return M
def word_fn(ws):
    fs = [T(w) for w in ws]
    def g(v):
        for f in reversed(fs): v = f(v)
        return v
    return g
def x_(i): return 1 << (i - 1)
def z_(i): return 1 << (2 + i)
def ee(i): return x_(i) | z_(i)
