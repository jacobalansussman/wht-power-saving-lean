# Independent recomputation (mod p, exact arithmetic) of ex-dimension's block-2 Jacobian ranks of the program map
#   u_j = rho(W_1) sum_z M_{j,z} z,  O_t = sum_z F_{t,z} z,  M, F in A = span rho(B),  inputs E_l = I_d in slot l.
# Mod-p rank at an F_p point of a sketch J S R (S random tangent directions, R random column projection) is a LOWER
# bound on the char-0 generic rank (all larger minors vanish identically over Z[i,1/2]_(p)).  Two primes.
import sys, os, json, time, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import letters, rep, rank_mod, Echelon, mm
def span_basis(mats, p):
    d = mats[0].shape[0]; E = Echelon(d * d, p)
    for c in range(0, len(mats), 256): E.add(np.stack([M.reshape(-1) for M in mats[c:c + 256]]))
    return E.rows.reshape(-1, d, d)
def group_elems(gens, p):
    n = gens[0].shape[0]; e = np.eye(n, dtype=np.int64); seen = {e.tobytes(): e}; fr = [e]
    while fr:
        nf = []
        for M in fr:
            for g in gens:
                X = mm(g, M, p); k = X.tobytes()
                if k not in seen: seen[k] = X; nf.append(X)
        fr = nf
    return list(seen.values())
def module_dim(gens_all, X0, p):
    """dim of span{w X0 : words w in gens_all} (closure), X0: (d, c)"""
    d, c = X0.shape; E = Echelon(d * c, p); fr = E.add(X0.reshape(1, -1))
    while fr:
        cand = np.concatenate([np.stack([mm(g, v.reshape(d, c), p).reshape(-1) for v in fr]) for g in gens_all])
        fr = E.add(cand)
    return E.dim()
class Prog:
    def __init__(self, m, repn, W, nsteps, p, seed):
        self.p = p; self.rng = np.random.default_rng(seed)
        free, Wr = letters(m, p)
        Bel = group_elems(free, p); self.nB = len(Bel)
        self.F = [rep(repn, g, p) for g in free]; self.W1 = rep(repn, Wr[1], p); self.K = rep(repn, Wr[m], p)
        self.A = span_basis([rep(repn, g, p) for g in Bel], p); self.s = self.A.shape[0]
        self.d = d = self.W1.shape[0]
        self.gU = module_dim(self.F + [self.W1], np.eye(d, dtype=np.int64), p)
        self.W = W; self.n = nsteps; self.slots = [W + j for j in range(nsteps)] + [W + nsteps] * W
        self.P = sum(self.slots) * self.s; self.Eeff = W * W * self.gU
    def inputs(self):
        d, W = self.d, self.W; out = []
        for l in range(W):
            X = np.zeros((d, d * W), np.int64); X[:, l * d:(l + 1) * d] = np.eye(d, dtype=np.int64); out.append(X)
        return out
    def rand_A(self, k):
        c = self.rng.integers(0, self.p, size=(k, self.s))
        return (c @ self.A.reshape(self.s, -1) % self.p).reshape(k, self.d, self.d)
    def random_point(self): return [list(self.rand_A(ns)) for ns in self.slots]
    def forward(self, Ms):
        p = self.p; src = self.inputs(); k = 0
        for j in range(self.n):
            c = sum(mm(Ms[k][i], src[i], p) for i in range(len(src))) % p; k += 1
            src.append(mm(self.W1, c, p))
        outs = []
        for t in range(self.W):
            outs.append(sum(mm(Ms[k][i], src[i], p) for i in range(len(src))) % p); k += 1
        return src, outs
    def jac_rank(self, Ms, extra=24):
        p = self.p; src, _ = self.forward(Ms); kdir = self.Eeff + extra
        d, c = self.d, self.d * self.W
        dsrc = [np.zeros((kdir, d, c), np.int64) for _ in range(self.W)]; k = 0
        for j in range(self.n):
            dc = np.zeros((kdir, d, c), np.int64)
            for i in range(len(dsrc)):
                D = self.rand_A(kdir)
                dc = (dc + np.matmul(D, src[i]) % p + np.matmul(Ms[k][i], dsrc[i]) % p) % p
            k += 1; dsrc.append(np.matmul(self.W1, dc) % p)
        outs = []
        for t in range(self.W):
            do = np.zeros((kdir, d, c), np.int64)
            for i in range(len(dsrc)):
                D = self.rand_A(kdir)
                do = (do + np.matmul(D, src[i]) % p + np.matmul(Ms[k][i], dsrc[i]) % p) % p
            k += 1; outs.append(do.reshape(kdir, -1))
        J = np.concatenate(outs, 1)                      # kdir x (W d c)
        R = self.rng.integers(0, p, size=(J.shape[1], kdir))
        JR = np.zeros((kdir, kdir), np.int64)
        for a in range(0, J.shape[1], 512): JR = (JR + J[:, a:a + 512] @ R[a:a + 512] % p) % p
        return rank_mod(JR, p)
if __name__ == '__main__':
    m, repn, W, n, p, seed = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    t0 = time.time(); Pr = Prog(m, repn, W, n, p, seed)
    r = Pr.jac_rank(Pr.random_point())
    print(json.dumps({'m': m, 'rep': repn, 'shape': f'{W}in{n}', 'p': p, '|B|': Pr.nB, 'dimA': Pr.s, 'gU': Pr.gU,
                      'P': Pr.P, 'Eeff': Pr.Eeff, 'rank_modp': r, 'codim_lowerbound_is': Pr.Eeff - r, 'sec': round(time.time() - t0, 1)}), flush=True)
