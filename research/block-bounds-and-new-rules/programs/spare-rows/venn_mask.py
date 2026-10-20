# venn_mask.py (d-G3-SPECTATOR/more-spectators, wave 17).  Two-factor Venn certificates W = Lambda^a V (x) Lambda^b V,
# V = C^(2^m), N_J = span{e_T1 (x) e_T2 : |T1 cap T2| in J}, for SPECTATOR targets B_k (x) I_(m-k) = dir(e1)...dir(ek).
#   * f(target) by the identity of wave16/b4-cert/demcalc.py (copied, generalised from the kernel to a bit mask):
#     g = B_k (x) I has g^2 = shift(e1+..+ek) (free, keeps N), Walsh eigenvalue i^{wt(s & mask)}; so
#     f(g) = sum_c rank(C[rows not in class c]) - 3 rank(C), C = Walsh coordinates of N (integer minors of the Hadamard).
#     Blocks: X-characters alpha (classes under bit permutations that keep the mask), and for even f the Z-shifts s with
#     f*wt(s&mask) + 2*wt(s&mask&alpha) = 0 mod 4.
#   * ranks: 'float' (SVD screen) or 'modp' (blocked elimination over GF(p), p = 999983 < 2^20, float64 BLAS exact):
#     a mod-p rank is <= the rational rank, and rank(C) = number of columns is asserted, so the mod-p value of
#     f(target) is a LOWER bound on the true f(target)  [the right direction for a lower bound on moves].
#   * cap = f(dir e1) EXACT (integer cut, wave16/b4-cert/venn2cut.py: dir(e1) is monomial in the pair basis).
# Bound on moves per array = f(target) / cap   (Theorem R; all rank-1 moves are conjugate to dir(e1) by free address maps).
# Usage: python3 -I venn_mask.py M A B K J RANK      e.g.  venn_mask.py 4 2 3 2 1 modp     (J like 1 or 0,2)
import sys, itertools, time
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import numpy as np
from math import comb, factorial

P = 999983
def popc(x): return bin(x).count('1')

# ---------------- exact cap (copied logic of venn2cut.py) ----------------
def cap_matrix(m, a, b):
    n = 2 ** m; Pn = n // 2
    states = [(i, j) for i in range(3) for j in range(3)]
    def multinom(ks):
        r = factorial(sum(ks))
        for k in ks: r //= factorial(k)
        return r
    def comps(total, parts):
        if parts == 1: yield (total,); return
        for k in range(total + 1):
            for rest in comps(total - k, parts - 1): yield (k,) + rest
    M = {}
    for cnt in comps(Pn, len(states)):
        if sum(c * s[0] for c, s in zip(cnt, states)) != a: continue
        if sum(c * s[1] for c, s in zip(cnt, states)) != b: continue
        d = dict(zip(states, cnt)); mult = multinom(cnt)
        other = 2 ** (d[(1, 0)] + d[(0, 1)] + d[(2, 1)] + d[(1, 2)])
        c0 = d[(2, 1)] + d[(1, 2)] + 2 * d[(2, 2)]; j = d[(1, 1)]
        for dp, dm, op in itertools.product(range(j + 1), repeat=3):
            om = j - dp - dm - op
            if om < 0: continue
            w = mult * other * multinom((dp, dm, op, om))
            c = c0 + dp + dm; c2 = c0 + op + dm
            M[(c, c2)] = M.get((c, c2), 0) + w
    return M
def exact_cap(m, a, b, J):
    M = cap_matrix(m, a, b)
    return sum(v for (x, y), v in M.items() if x in J and y not in J)

# ---------------- fast rank mod p ----------------
def rank_modp_fast(A, p=P, panel=64):
    A = np.asarray(A, dtype=np.int64) % p
    rank = 0
    while A.shape[0] > 0 and A.shape[1] > 0:
        R, Cn = A.shape; b = min(panel, Cn)
        Pn = A[:, :b].copy(); rowid = np.arange(R); r = 0; pr = []; pc = []
        for j in range(b):
            if r == R: break
            nz = np.nonzero(Pn[r:, j])[0]
            if len(nz) == 0: continue
            i = r + nz[0]
            if i != r:
                Pn[[r, i]] = Pn[[i, r]]; rowid[[r, i]] = rowid[[i, r]]
            Pn[r] = (Pn[r] * pow(int(Pn[r, j]), p - 2, p)) % p
            col = Pn[r + 1:, j].copy()
            nzr = np.nonzero(col)[0]
            if len(nzr):
                Pn[r + 1 + nzr] = (Pn[r + 1 + nzr] - (col[nzr, None] * Pn[r][None, :]) % p) % p
            pr.append(int(rowid[r])); pc.append(j); r += 1
        t = len(pr); rank += t
        if t == 0:
            A = A[:, b:]; continue
        others = np.setdiff1d(np.arange(R), np.array(pr))
        if len(others) == 0 or Cn == b: break
        Pv = A[np.ix_(pr, pc)]
        Pinv = inv_modp(Pv, p)
        X = A[np.ix_(pr, np.arange(b, Cn))]
        Mx = matmul_modp(Pinv, X, p)
        Y = A[np.ix_(others, pc)]
        Z = A[np.ix_(others, np.arange(b, Cn))]
        A = (Z - matmul_modp(Y, Mx, p)) % p
    return rank
def matmul_modp(X, Y, p=P):
    # entries < 2^20, inner dimension chunks of 64 -> partial sums < 2^46, exact in float64
    X = np.asarray(X, dtype=np.float64); Y = np.asarray(Y, dtype=np.float64)
    out = np.zeros((X.shape[0], Y.shape[1]), dtype=np.float64)
    for s in range(0, X.shape[1], 64):
        out = np.fmod(out + X[:, s:s + 64] @ Y[s:s + 64, :], p)
    return out.astype(np.int64) % p
def inv_modp(Mt, p=P):
    t = Mt.shape[0]; Aug = np.hstack([Mt % p, np.eye(t, dtype=np.int64)])
    for c in range(t):
        piv = c + np.nonzero(Aug[c:, c])[0][0]
        if piv != c: Aug[[c, piv]] = Aug[[piv, c]]
        Aug[c] = (Aug[c] * pow(int(Aug[c, c]), p - 2, p)) % p
        col = Aug[:, c].copy(); col[c] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): Aug[nzr] = (Aug[nzr] - (col[nzr, None] * Aug[c][None, :]) % p) % p
    return Aug[:, t:]

# ---------------- f(target) for a mask ----------------
class Demand:
    def __init__(self, m, avec, inS, mask, rankmode='float'):
        self.m, self.avec, self.r = m, tuple(avec), len(avec)
        self.n = n = 2 ** m; self.f = sum(avec); self.mask = mask; self.rankmode = rankmode
        self.subs = {a: list(itertools.combinations(range(n), a)) for a in set(avec)}
        self.sidx = {a: {T: i for i, T in enumerate(self.subs[a])} for a in set(avec)}
        Hm = np.array([[(-1) ** popc(s & x) for x in range(n)] for s in range(n)], float)
        self.lab = [popc(s & mask) for s in range(n)]
        self.Hm = Hm; self.comp = {}
        for a in set(avec):
            S = self.subs[a]; idxU = np.array(S)
            if len(S) ** 2 > 2e7:          # too big to store (> 160 MB): minors computed on the fly (demcalc.py's method)
                self.comp[a] = None; continue
            M = np.zeros((len(S), len(S)))
            for i, U in enumerate(S):
                sub = Hm[np.array(U)][:, idxU]; sub = np.transpose(sub, (1, 0, 2))
                M[i] = np.linalg.det(sub)
            self.comp[a] = np.rint(M).astype(np.int64)
        self.inS = inS
    def minors(self, a, ri, ci):
        Sa = np.array(self.subs[a]); out = np.empty((len(ri), len(ci)), dtype=np.int64)
        Tc = Sa[ci]
        ch = max(1, int(4e6 // max(1, len(ci) * a * a)))
        for st in range(0, len(ri), ch):
            Ur = Sa[ri[st:st + ch]]
            blk = self.Hm[Ur[:, None, :, None], Tc[None, :, None, :]]
            out[st:st + ch] = np.rint(np.linalg.det(blk)).astype(np.int64)
        return out
    def sigma(self, tup):
        s = 0
        for T in tup:
            for x in T: s ^= x
        return s
    def run(self, verbose=False):
        m, n, f, mask = self.m, self.n, self.f, self.mask
        allt = list(itertools.product(*[self.subs[a] for a in self.avec]))
        St = [t for t in allt if self.inS(t)]
        dimN = len(St)
        def shift(t, z): return tuple(tuple(sorted(x ^ z for x in T)) for T in t)
        St = sorted(set(min(shift(t, z) for z in range(n)) for t in St))
        tsig = np.array([self.sigma(t) for t in St])
        tidx = [np.array([self.sidx[a][t[i]] for t in St]) for i, a in enumerate(self.avec)]
        usig = np.array([self.sigma(u) for u in allt])
        ucls = np.array([sum(self.lab[s] for S in u for s in S) % 4 for u in allt])
        uidx = [np.array([self.sidx[a][u[i]] for u in allt]) for i, a in enumerate(self.avec)]
        inb = [i for i in range(m) if (mask >> i) & 1]; outb = [i for i in range(m) if not (mask >> i) & 1]
        alphas = []
        for w1 in range(len(inb) + 1):
            for w2 in range(len(outb) + 1):
                al = sum(1 << i for i in inb[:w1]) | sum(1 << i for i in outb[:w2])
                alphas.append((comb(len(inb), w1) * comb(len(outb), w2), al))
        pos = {u: i for i, u in enumerate(allt)}
        total = 0; dimchk = 0
        for (wt_a, alpha) in alphas:
            rows0 = np.nonzero(usig == alpha)[0]
            Z0 = [s for s in range(n) if (f * popc(s & mask) + 2 * popc(alpha & s & mask)) % 4 == 0] if f % 2 == 0 else [0]
            def key_of(sig): return tuple(popc(s & sig) & 1 for s in Z0)
            rset = set()
            for i in rows0:
                u = allt[i]
                rset.add(min(pos[tuple(tuple(sorted(x ^ s_ for x in S)) for S in u)] for s_ in Z0))
            rows = np.array(sorted(rset))
            colkeys = {}
            for j, sg in enumerate(tsig): colkeys.setdefault(key_of(int(sg)), []).append(j)
            blk_tot = 0
            for kk, cols in colkeys.items():
                cols = np.array(cols)
                C = np.ones((len(rows), len(cols)), dtype=np.int64)
                for i in range(self.r):
                    a = self.avec[i]
                    if self.comp[a] is not None:
                        C *= self.comp[a][np.ix_(uidx[i][rows], tidx[i][cols])]
                    else:
                        C *= self.minors(a, uidx[i][rows], tidx[i][cols])
                rk = self.rank(C)      # orbit representatives may project to 0 in this character: rk <= len(cols)
                rc = [self.rank(C[ucls[rows] != c]) for c in range(4)]
                blk = sum(rc) - 3 * rk
                blk_tot += blk; dimchk += wt_a * rk
            total += wt_a * blk_tot
            if verbose: print(f"  alpha={alpha:0{m}b} x{wt_a}: rows {len(rows)} part {blk_tot}", flush=True)
        # sum over blocks of rank(C) = dim N exactly: with mod-p ranks (each <= the rational rank) this forces every
        # rank(C) to be exact, so the mod-p total is a LOWER bound on f(target) (the rank(C[rows not in c]) terms can only drop)
        assert dimchk == dimN, (dimchk, dimN)
        return dimN, total
    def rank(self, C):
        if C.shape[0] == 0 or C.shape[1] == 0: return 0
        if self.rankmode == 'modp':
            # put the short side as columns
            return rank_modp_fast(C if C.shape[0] >= C.shape[1] else C.T)
        s = np.linalg.svd(C.astype(float), compute_uv=False)
        return int((s > 1e-7 * max(1.0, s[0])).sum())

def venn_pred(J):
    def pred(t):
        return len(set(t[0]) & set(t[1])) in J
    return pred

if __name__ == '__main__':
    m, a, b, k = (int(x) for x in sys.argv[1:5]); J = [int(x) for x in sys.argv[5].split(',')]
    rk = sys.argv[6] if len(sys.argv) > 6 else 'float'
    t0 = time.time()
    cap = exact_cap(m, a, b, J)
    mask = (1 << k) - 1
    D = Demand(m, [a, b], venn_pred(J), mask, rankmode=rk)
    dN, val = D.run(verbose=True)
    print(f"m={m} Lambda^{a} (x) Lambda^{b} J={J}: dim N {dN}; cap (exact) {cap}; f(B_{k} (x) I_{m-k}) [{rk}] = {val}; "
          f"bound {val}/{cap} = {val/cap:.6f} moves per array  ({time.time()-t0:.1f}s)", flush=True)
