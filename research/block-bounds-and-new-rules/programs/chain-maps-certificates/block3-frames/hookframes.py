"""xc3-cert3: the block-3 hook certificate is a PAULI-BLIND frame family on the 135 rational frames.

Y_U = (g (x) Lambda^4 g) N  for any Clifford g with g L0 = U, where N = span{e_x (x) e_T : x in T, |T| = 4} in
C^8 (x) Lambda^4 C^8 (d-tof-cert's hook, dim 280).  Well defined: Stab(L0) in the Clifford group = monomial Cliffords
(Paulis included), and every monomial m keeps N.  So delta(U, V) = dim(Y_U + Y_V) - 280 depends only on the Clifford
orbit of (U, V), i.e. on the distance, and equals f(g) of the hook certificate: 80 (direction), 144 (rank-2 block),
224 (kernel).  The family is indexed by frames, so with the Averaging Lemma it is a Pauli-blind certificate
(Paulis fix every frame).  Check here [floats]: BFS over Clifford generators (H, S on each qubit, CNOTs) from L0,
keeping one Clifford unitary per frame; 135 frames; delta(L0, V) for all V, and delta(U, V) for random pairs U, V.
"""
import itertools, sys
import numpy as np

n = 8
def kron(*ms):
    out = np.ones((1, 1))
    for m in ms: out = np.kron(out, m)
    return out
I2 = np.eye(2); Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2); Sg = np.diag([1, 1j])
def on(q, m):
    ms = [I2] * 3; ms[q] = m; return kron(*ms)
def cnot(c, t):
    M = np.zeros((8, 8))
    for x in range(8):
        bits = [(x >> (2 - q)) & 1 for q in range(3)]
        if bits[c]: bits[t] ^= 1
        y = sum(b << (2 - q) for q, b in enumerate(bits)); M[y, x] = 1
    return M
GENS = [on(q, Hd) for q in range(3)] + [on(q, Sg) for q in range(3)] + [cnot(c, t) for c in range(3) for t in range(3) if c != t]

# Pauli labels: P(x|z) = X^x Z^z ; frame of a Clifford g = { u : g P_z g^-1 is diagonal-type } = g L0 (labels)
def pauli(x, z):
    ms = []
    for q in range(3):
        m = np.eye(2)
        if (x >> (2 - q)) & 1: m = np.array([[0, 1], [1, 0]]) @ m
        if (z >> (2 - q)) & 1: m = m @ np.diag([1, -1])
        ms.append(m)
    return kron(*ms)
PAULIS = {(x, z): pauli(x, z) for x in range(8) for z in range(8)}
def frame_of(g):
    lab = []
    for z in range(1, 8):
        M = g @ PAULIS[(0, z)] @ g.conj().T
        for (x2, z2), Pm in PAULIS.items():
            c = np.trace(Pm.conj().T @ M) / 8
            if abs(abs(c) - 1) < 1e-9: lab.append((x2, z2)); break
    return frozenset(lab)

subs4 = list(itertools.combinations(range(8), 4)); s4 = {T: i for i, T in enumerate(subs4)}
Ncols = [(x, s4[T]) for T in subs4 for x in T]
def compound(g, k):
    subs = list(itertools.combinations(range(8), k))
    C = np.zeros((len(subs), len(subs)), complex)
    for i, T in enumerate(subs):
        for j, S in enumerate(subs):
            C[i, j] = np.linalg.det(g[np.ix_(T, S)])
    return C
def Ymat(g):
    C4 = compound(g, 4)
    big = np.kron(g, C4)                   # index (x, T) -> x * 70 + T
    cols = [x * 70 + t for x, t in Ncols]
    return big[:, cols]
def rank(M, tol=1e-8):
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > tol * s[0]).sum()), (s[s > tol * s[0]].min() / s[0], (s[s <= tol * s[0]].max() / s[0]) if (s <= tol * s[0]).any() else 0)

if __name__ == "__main__":
    g0 = np.eye(8, dtype=complex)
    L0 = frame_of(g0)
    rep = {L0: g0}; fr = [L0]
    while fr:
        nxt = []
        for U in fr:
            for G in GENS:
                g = G @ rep[U]; V = frame_of(g)
                if V not in rep: rep[V] = g; nxt.append(V)
        fr = nxt
    frames = list(rep)
    print("frames:", len(frames), flush=True)
    def dist(U, V):
        # labels as sets of 7 nonzero labels; |U cap V| = 2^(3-d) - 1
        return {7: 0, 3: 1, 1: 2, 0: 3}[len(U & V)]
    Y = {U: Ymat(rep[U]) for U in frames}
    print("dim Y_L0 =", rank(Y[L0])[0], flush=True)
    from collections import defaultdict, Counter
    res = defaultdict(Counter); worst = 1.0
    for V in frames:
        r, gap = rank(np.hstack([Y[L0], Y[V]])); res[dist(L0, V)][r - 280] += 1; worst = min(worst, gap[0] / max(gap[1], 1e-300))
    print("delta(L0, V) by distance:", {d: dict(c) for d, c in sorted(res.items())})
    rng = np.random.default_rng(1)
    res2 = defaultdict(Counter)
    for _ in range(int(sys.argv[1]) if len(sys.argv) > 1 else 60):
        U, V = frames[rng.integers(135)], frames[rng.integers(135)]
        r, gap = rank(np.hstack([Y[U], Y[V]])); res2[dist(U, V)][r - 280] += 1
    print("delta(U, V) on random pairs by distance:", {d: dict(c) for d, c in sorted(res2.items())})
    print("ratio delta(distance 3)/delta(distance 1) =", res[3].most_common(1)[0][0], "/", res[1].most_common(1)[0][0])
