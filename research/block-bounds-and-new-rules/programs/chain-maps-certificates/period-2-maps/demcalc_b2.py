# b4-cert (wave 16).  EXACT-ENTRY demand (and cap) of a coordinate certificate N = span{e_t : t in S} in
# W = Lambda^{a_1} V (x) ... (x) Lambda^{a_r} V, V = C^(2^m).
# Facts used: K^2 = shift(1...1) is free and keeps N, so N cap K N is K-invariant, i.e. N cap KN = (+)_c (N cap W_c),
# W_c = K-eigenspace i^c = span of Walsh tensors h_u (h_s = sum_x (-1)^{s.x} e_x; u = (S_1..S_r)) with sum of wt(s) = c
# mod 4.  Hence dem = sum_c rank(C[rows not in c]) - 3 rank(C), C[u,t] = <h_u, e_t> = prod_i det[(-1)^{s.x}]_{S_i,T_i}
# (integers).  Same for any element g with g^2 free: cap uses dir(e1) (eigenbasis: Walsh in bit 1 only).
# Blocks: X-shift characters alpha (K commutes with X_z; h_u has character sigma(u) = xor of all labels), and for
# even f the Z-shifts s with i^{f wt(s)} (-1)^{alpha.s} = 1 (they keep every W_c inside E_alpha); bit permutations
# (commute with K, keep N) identify alpha of equal weight.
# Ranks: float SVD (screening) or mod p (exact lower bound for each rank).
import sys, itertools, time
# (needs numpy; the original loaded a private copy of numpy and scipy here)
# xc-algebra copy of wave16/b4-cert/demcalc.py with a B2 (rank-2 block) mode added; nothing else changed.
import numpy as np
from math import comb

def popc(x): return bin(x).count('1')
_COMP = {}

class Demand:
    def __init__(self, m, avec, inS, modp=None, which='K'):
        """inS(t) -> bool for a standard tuple t = (T_1..T_r) of sorted tuples"""
        self.m, self.avec, self.r = m, tuple(avec), len(avec)
        self.n = n = 2 ** m; self.f = sum(avec); self.modp = modp; self.which = which
        self.subs = {a: list(itertools.combinations(range(n), a)) for a in set(avec)}
        self.sidx = {a: {T: i for i, T in enumerate(self.subs[a])} for a in set(avec)}
        # compound matrices of the (unnormalised) Hadamard  H[s,x] = (-1)^{s.x}   (for 'K')
        # for 'R' (dir e1): eigenbasis = Walsh in bit 0 only: g_s = e_s + (-1)^{s_0} e_{s^1}  (for s with...): use
        # basis vectors w_y = e_{y} + e_{y^1} (y even) and e_{y^1}... -> build matrix E[s, x] below.
        Hm = np.array([[(-1) ** popc(s & x) for x in range(n)] for s in range(n)], float)
        # K = dir(e1)...dir(em): Walsh h_s has eigenvalue i^{wt s};  R = dir(e1): eigenvalue i^{s_0} (bit 0 of s)
        # xc-algebra: 'B2' = rank-2 block dir(e1)dir(e2) on bits 0,1: Walsh h_s has eigenvalue i^{s0+s1}
        self.mask = {'K': n - 1, 'R': 1, 'B2': 3}[which]
        self.lab = [popc(s & self.mask) for s in range(n)]
        self.Hm = Hm
        self.comp = {}
        global _COMP
        for a in set(avec):
            if (m, a) in _COMP:
                self.comp[a] = _COMP[(m, a)]; continue
            S = self.subs[a]
            if comb(n, a) ** 2 <= 4e7:
                M = np.array([[np.linalg.det(Hm[np.ix_(U, T)]) if a else 1.0 for T in S] for U in S]) if a <= 2 else None
                if M is None:
                    idxU = np.array(S); 
                    M = np.zeros((len(S), len(S)))
                    for i, U in enumerate(S):
                        sub = Hm[np.array(U)][:, idxU]          # (a, nT, a)
                        sub = np.transpose(sub, (1, 0, 2))
                        M[i] = np.linalg.det(sub)
                self.comp[a] = np.rint(M).astype(np.int64); _COMP[(m, a)] = self.comp[a]
            else:
                self.comp[a] = None          # computed on the fly (chunked batched determinants)
        self.inS = inS

    def sigma(self, tup):
        s = 0
        for T in tup:
            for x in T: s ^= x
        return s

    def run(self, verbose=False):
        m, n, f = self.m, self.n, self.f
        allt = list(itertools.product(*[self.subs[a] for a in self.avec]))
        St = [t for t in allt if self.inS(t)]
        dimN = len(St)
        def shift(t, z): return tuple(tuple(sorted(x ^ z for x in T)) for T in t)
        St = sorted(set(min(shift(t, z) for z in range(n)) for t in St))
        tsig = np.array([self.sigma(t) for t in St])
        tidx = [np.array([self.sidx[a][T] for t in St for T in [t[i]]]) for i, a in enumerate(self.avec)]
        # Walsh tuples: same index set allt; label sigma and class
        usig = np.array([self.sigma(u) for u in allt])
        ucls = np.array([sum(self.lab[s] for S in u for s in S) % 4 for u in allt])
        uidx = [np.array([self.sidx[a][u[i]] for u in allt]) for i, a in enumerate(self.avec)]
        if self.which == 'K':
            alphas = [(1, (1 << w) - 1) for w in range(m + 1)]
            alphas = [(comb(m, w), (1 << w) - 1) for w in range(m + 1)]
        elif self.which == 'R':
            # bit permutations fixing bit 0 commute with R: alpha classes by (alpha_0, weight of the other bits)
            alphas = [(comb(m - 1, w), b0 | (((1 << w) - 1) << 1)) for b0 in (0, 1) for w in range(m)]
        else:
            # B2: bit permutations keeping {0,1} setwise commute with dir(e1)dir(e2)
            alphas = [(comb(2, w1) * comb(m - 2, w2), ((1 << w1) - 1) | (((1 << w2) - 1) << 2))
                      for w1 in range(3) for w2 in range(m - 1)]
        total = 0; dimchk = 0
        only = getattr(self, 'only', None)          # xc-algebra: run a single alpha class (for parallel runs)
        for ai, (wt_a, alpha) in enumerate(alphas):
            if only is not None and ai != only: continue
            rows0 = np.nonzero(usig == alpha)[0]
            # Z subgroup keeping the classes: s with (f*wt(s) + 2 alpha.s) % 4 == 0  (K);  for R: f*s0 + 2 alpha.s
            if self.which == 'K':
                Z0 = [s for s in range(n) if (f * popc(s) + 2 * popc(alpha & s)) % 4 == 0] if f % 2 == 0 else [0]
            else:
                mk = self.mask
                Z0 = [s for s in range(n) if (f * popc(s & mk) + 2 * popc(alpha & s & mk)) % 4 == 0] if f % 2 == 0 else [0]
            # characters of Z0 on standard side: chi_t(s) = (-1)^{s.sigma(t)}  -> class key = tuple of parities
            def key_of(sig): return tuple(popc(s & sig) & 1 for s in Z0)
            pos = {u: i for i, u in enumerate(allt)} if not hasattr(self, '_pos') else self._pos
            self._pos = pos
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
                # rows: all Walsh tuples with sigma = alpha (no orbit reduction: rank unchanged, cost higher)
                C = np.ones((len(rows), len(cols)), dtype=np.int64)
                for i in range(self.r):
                    a = self.avec[i]
                    if self.comp[a] is not None:
                        C *= self.comp[a][np.ix_(uidx[i][rows], tidx[i][cols])]
                    else:
                        C *= self.minors(a, uidx[i][rows], tidx[i][cols])
                rk = self.rank(C)
                rc = [self.rank(C[ucls[rows] != c]) for c in range(4)]
                blk = sum(rc) - 3 * rk
                blk_tot += blk; dimchk += wt_a * rk
            total += wt_a * blk_tot
            if verbose: print(f"  alpha={alpha} x{wt_a}: part {blk_tot}", flush=True)
        if getattr(self, 'only', None) is None:
            assert dimchk == dimN, (dimchk, dimN)   # orbit reduction must not lose dimension
        else:
            print(f"  class {self.only}: weighted rank {dimchk} (sum over classes must be dim N {dimN})", flush=True)
        return dimN, total

    def minors(self, a, ri, ci):
        Sa = np.array(self.subs[a]); out = np.empty((len(ri), len(ci)), dtype=np.int64)
        Hm = self.Hm
        Tc = Sa[ci]                                   # (C, a)
        ch = max(1, int(2e7 // max(1, len(ci) * a * a)))
        for st in range(0, len(ri), ch):
            Ur = Sa[ri[st:st + ch]]                   # (R', a)
            blk = Hm[Ur[:, None, :, None], Tc[None, :, None, :]]   # (R', C, a, a)
            out[st:st + ch] = np.rint(np.linalg.det(blk)).astype(np.int64)
        return out

    def rank(self, C):
        if C.shape[0] == 0 or C.shape[1] == 0: return 0
        if self.modp:
            return rank_modp(C % self.modp, self.modp)
        s = np.linalg.svd(C.astype(float), compute_uv=False)
        return int((s > 1e-7 * max(1.0, s[0])).sum())

def rank_modp(M, p):
    M = M.copy() % p; r = 0; rows, cols = M.shape
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        piv = r + nz[0]
        if piv != r: M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): M[nzr] = (M[nzr] - (col[nzr, None] * M[r][None, :]) % p) % p
        r += 1
    return r

def venn_vec(t, r):
    cnt = {}
    for x in set(x for T in t for x in T):
        reg = tuple(i for i in range(r) if x in t[i])
        cnt[reg] = cnt.get(reg, 0) + 1
    return cnt

if __name__ == '__main__':
    m = int(sys.argv[1]); avec = [int(x) for x in sys.argv[2].split(',')]
    which = sys.argv[3]            # 'K', 'R' or 'B2'
    pred = sys.argv[4]             # python expression in I (dict: tuple(region) -> count), e.g. "I.get((0,1),0)==1"
    t0 = time.time()
    r = len(avec)
    D = Demand(m, avec, lambda t: eval(pred, {}, {'I': venn_vec(t, r)}), which=which)
    if len(sys.argv) > 5: D.only = int(sys.argv[5])
    dN, val = D.run(verbose=True)
    print(f"m={m} {avec} N: {pred}: dim N {dN}, f({which}) = {val}  ({time.time()-t0:.1f}s)")
