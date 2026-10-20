"""check_orth.py -- v-w7-xwborder (skeptic).  OWN CODE (no xw2lib).  Exact normalised traces of the orthogonality facts
Theorem 2W / 2W' rest on, with a free menu WIDER than lemD.py's: column-uniform Paulis, CNOTs, S, CS, CZ; chain maps
(x_c ^= x_{c-1}) and their prefix-XOR inverses; Toffoli across columns; parity phases Pr_r = a + b Zs_r; and
NON-monomial bounded-support boundary letters (one-column H, one-column sqrt-X = dir, one-column CS·H).

Facts tested (tau = tr/N, computed exactly as the sum of all N diagonal entries, by applying the operator to basis
vectors in batches):
 (O')  block 2: |tau(w^* C D_z C^-1)| for free words w, C          -> must decay (bound (s c)^L 2^(-f/2))
 (O1)  block 2: |tau((c K c')^* c1 D_z c2)|  (Lemma O, j = 1)       -> must decay like 2^(-f/2)
 (R1)  block 3: |tau((c R c')^* c1 D_z c2)|, R = (dir e1 dir e2)^(x)f (2W' for rank-2 targets) -> must decay
 (K2)  block 3: |tau((c K c')^* c1 D_z c2 D_y c3)|  (Lemma O, j = 2) -> must decay
 controls: tau(K^* D1 D2) = 1 at block 2; tau(R^* D1 D2) = 1 at block 3 (no decay at the top layer).
usage: python3 check_orth.py m fmax samples seed
"""
import sys, itertools
import numpy as np
from functools import reduce
def prodb(xs):
    return reduce(lambda p, q: p * q, xs)

A = (1 + 1j) / 2; B = (1 - 1j) / 2

class Game:
    def __init__(self, m, f):
        self.m, self.f, self.nb = m, f, m * f
        self.N = 2 ** self.nb
    # vectors: shape (N, batch); bit (c, r) is axis c*m + r of the reshaped tensor
    def T(self, V):
        return V.reshape([2] * self.nb + [V.shape[1]])
    def U(self, T):
        return T.reshape(self.N, -1)
    def ax(self, c, r):
        return c * self.m + r
    # ---- elementary operations on tensors ----
    def flip(self, T, bits):
        for (c, r) in bits:
            T = np.flip(T, axis=self.ax(c, r))
        return T
    def bitval(self, c, r):
        shp = [1] * (self.nb + 1); shp[self.ax(c, r)] = 2
        return np.arange(2).reshape(shp)
    def phase(self, T, fn):        # diagonal: multiply by fn(bit-value arrays)
        return T * fn(self)
    def perm(self, V, mapfn):      # monomial permutation: new[idx'] = old[idx], idx' = mapfn(bits of idx)
        N = self.N; idx = np.arange(N)
        bits = [(idx >> (self.nb - 1 - k)) & 1 for k in range(self.nb)]   # axis k <-> bit nb-1-k of the flat index
        nb = mapfn([b.copy() for b in bits])
        new = np.zeros_like(idx)
        for k in range(self.nb):
            new |= nb[k] << (self.nb - 1 - k)
        out = np.empty_like(V); out[new] = V
        return out

def dirz_col(G, T, c, z):          # apply dir(z) = a I + b X^z on column c
    bits = [(c, r) for r in range(G.m) if (z >> r) & 1]
    return A * T + B * G.flip(T, bits)

def D(G, V, z, inv=False):         # paid letter dir(z)^(x)f (inverse: X^z dir(z))
    T = G.T(V)
    for c in range(G.f):
        T = dirz_col(G, T, c, z)
    if inv:
        T = G.flip(T, [(c, r) for c in range(G.f) for r in range(G.m) if (z >> r) & 1])
    return G.U(T)

# ---- free letters: name -> (apply, apply_adjoint); all bounded support uniformly in f ----
def letters(G):
    m, f = G.m, G.f
    L = {}
    def diag(fn):
        return (lambda V: G.U(G.T(V) * fn()), lambda V: G.U(G.T(V) * np.conj(fn())))
    for r in range(m):
        L["X%d" % r] = (lambda V, r=r: G.U(G.flip(G.T(V), [(c, r) for c in range(f)])),) * 2
        L["Z%d" % r] = diag(lambda r=r: prodb([(-1.0) ** G.bitval(c, r) for c in range(f)]))
        L["S%d" % r] = diag(lambda r=r: prodb([(1j) ** G.bitval(c, r) for c in range(f)]))
        # parity phase Pr_r = a + b Zs_r
        L["Pr%d" % r] = diag(lambda r=r: A + B * prodb([(-1.0) ** G.bitval(c, r) for c in range(f)]))
    if m >= 2:
        L["CS"] = diag(lambda: prodb([(1j) ** (G.bitval(c, 0) * G.bitval(c, 1)) for c in range(f)]))
        L["CZ"] = diag(lambda: prodb([(-1.0) ** (G.bitval(c, 0) * G.bitval(c, 1)) for c in range(f)]))
        def cnot(bits):
            for c in range(f):
                bits[G.ax(c, 1)] = bits[G.ax(c, 1)] ^ bits[G.ax(c, 0)]
            return bits
        def cnot_inv(bits):
            return cnot(bits)
        L["C01"] = (lambda V: G.perm(V, cnot), lambda V: G.perm(V, cnot_inv))
        # Toffoli across columns: (c,0) ^= (c-1,0) & (c-1,1), old values (lower-triangular, invertible)
        def toff(bits):
            old = [b.copy() for b in bits]
            for c in range(1, f):
                bits[G.ax(c, 0)] = old[G.ax(c, 0)] ^ (old[G.ax(c - 1, 0)] & old[G.ax(c - 1, 1)])
            return bits
        def toff_inv(bits):
            out = [b.copy() for b in bits]
            for c in range(1, f):    # recover sequentially (x_{c-1} already recovered)
                out[G.ax(c, 0)] = bits[G.ax(c, 0)] ^ (out[G.ax(c - 1, 0)] & out[G.ax(c - 1, 1)])
            return out
        L["Tof"] = (lambda V: G.perm(V, toff), lambda V: G.perm(V, toff_inv))
    for r in range(min(m, 2)):
        def chain(bits, r=r):        # x_c ^= x_{c-1} (old values), row r
            old = [b.copy() for b in bits]
            for c in range(1, f):
                bits[G.ax(c, r)] = old[G.ax(c, r)] ^ old[G.ax(c - 1, r)]
            return bits
        def chain_inv(bits, r=r):    # prefix XOR
            out = [b.copy() for b in bits]
            for c in range(1, f):
                out[G.ax(c, r)] = bits[G.ax(c, r)] ^ out[G.ax(c - 1, r)]
            return out
        L["a%d" % r] = (lambda V, r=r: G.perm(V, lambda b: chain(b, r)), lambda V, r=r: G.perm(V, lambda b: chain_inv(b, r)))
        L["a%di" % r] = (L["a%d" % r][1], L["a%d" % r][0])
    # non-monomial boundary letters (one column only): H on (0,0); dir(e_{m-1}) on the last column; CS.H on column 0
    Hm = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    def hcol(V, adj=False):
        T = G.T(V); a = G.ax(0, 0)
        T = np.moveaxis(np.tensordot(Hm, T, axes=([1], [a])), 0, a)
        return G.U(T)
    L["fH"] = (hcol, hcol)
    zl = 1 << (m - 1)
    L["lDir"] = (lambda V: G.U(dirz_col(G, G.T(V), f - 1, zl)),
                 lambda V: G.U(np.conj(A) * G.T(V) + np.conj(B) * G.flip(G.T(V), [(f - 1, m - 1)])))
    return L

def apply_word(L, word, V, adj=False):
    # word = list of letter names applied right-to-left as matrices: w = L[w0] L[w1] ... ; w V = L[w0](L[w1](...V))
    if not adj:
        for nm in reversed(word):
            V = L[nm][0](V)
    else:                            # w^* V = L[w_last]^* ... L[w0]^* V
        for nm in word:
            V = L[nm][1](V)
    return V

def trace(G, op, batch=512):
    N = G.N; tot = 0j
    for s in range(0, N, batch):
        e = min(N, s + batch)
        E = np.zeros((N, e - s), complex); E[np.arange(s, e), np.arange(e - s)] = 1
        Y = op(E)
        tot += np.sum(Y[np.arange(s, e), np.arange(e - s)])
    return tot / N

def main():
    m, fmax, samples, seed = (int(t) for t in sys.argv[1:5])
    rng = np.random.default_rng(seed)
    for f in range(1, fmax + 1):
        G = Game(m, f); L = letters(G); names = sorted(L)
        def rw(maxlen):
            return [names[i] for i in rng.integers(0, len(names), rng.integers(0, maxlen + 1))]
        zs = list(range(1, 2 ** m))
        res = {}
        if m == 2:
            # (O') tau(w^* C D C^-1) : op(V) = w^* C D C^* V
            worst = 0
            special = [(["lDir"], [], 2), (["lDir"], [], 3), (["fH", "lDir"], [], 3), (["Pr0", "Pr1"], [], 1),
                       (["Pr1"], ["CS"], 3), ([], ["a0"], 1), (["Tof"], ["a1"], 3), (["lDir", "Pr1"], ["fH"], 2)]
            for t in range(samples + len(special)):
                if t < len(special): w, C, z = special[t]
                else: w, C, z = rw(3), rw(3), zs[rng.integers(len(zs))]
                op = lambda V: apply_word(L, w, apply_word(L, C, D(G, apply_word(L, C, V, adj=True), z)), adj=True)
                worst = max(worst, abs(trace(G, op)))
            res["O'"] = worst
            # (O1) tau((cKc')^* c1 D c2) = tau(c'^* K^* c^* c1 D c2)
            worst = 0
            for t in range(samples):
                c, cp, c1, c2, z = rw(2), rw(2), rw(2), rw(2), zs[rng.integers(len(zs))]
                inv = bool(rng.integers(2))
                def op(V):
                    V = apply_word(L, c2, V); V = D(G, V, z, inv); V = apply_word(L, c1, V)
                    V = apply_word(L, c, V, adj=True)
                    V = D(G, D(G, V, 2, inv=True), 1, inv=True)      # K^* = D2^* D1^* applied as (D1 D2)^*: K = D1 D2? K^*V
                    V = apply_word(L, cp, V, adj=True)
                    return V
                worst = max(worst, abs(trace(G, op)))
            res["O1"] = worst
            # control: tau(K^* D1 D2) with K = dir(e1) dir(e2) (commuting) -> 1
            ctl = trace(G, lambda V: D(G, D(G, D(G, D(G, V, 2), 1), 1, inv=True), 2, inv=True))
            res["ctl"] = abs(ctl)
        else:
            worst = 0; worstK = 0
            for t in range(samples):
                c, cp, c1, c2, z = rw(2), rw(2), rw(2), rw(2), zs[rng.integers(len(zs))]
                def op(V):
                    V = apply_word(L, c2, V); V = D(G, V, z); V = apply_word(L, c1, V)
                    V = apply_word(L, c, V, adj=True)
                    V = D(G, D(G, V, 2, inv=True), 1, inv=True)      # R^*
                    V = apply_word(L, cp, V, adj=True)
                    return V
                worst = max(worst, abs(trace(G, op)))
                y = zs[rng.integers(len(zs))]; c3 = rw(2)
                def opK(V):
                    V = apply_word(L, c3, V); V = D(G, V, y); V = apply_word(L, c2, V); V = D(G, V, z)
                    V = apply_word(L, c1, V); V = apply_word(L, c, V, adj=True)
                    for r in reversed(range(m)):
                        V = D(G, V, 1 << r, inv=True)
                    V = apply_word(L, cp, V, adj=True)
                    return V
                worstK = max(worstK, abs(trace(G, opK)))
            res["R1"] = worst; res["K2"] = worstK
            res["ctl"] = abs(trace(G, lambda V: D(G, D(G, D(G, D(G, V, 2), 1), 1, inv=True), 2, inv=True)))
        print("m=%d f=%d N=%d  " % (m, f, G.N) + "  ".join("%s %.4e" % kv for kv in res.items()) +
              "   2^(-f/2) %.4e" % 2 ** (-f / 2), flush=True)

if __name__ == "__main__":
    main()
