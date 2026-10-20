"""relax.py -- w5-blueprint: EXACT frame-sum relaxation of an H-symmetric block-4 skeleton with the top orbit's content
left completely free (every letter on every earlier generator), decided per character by linear algebra.

Model (sp8.py): a term [g] x_l is stored as (l, R) with R = sigma(g)^-1 (8 columns); its cell is R(L_Z).  A move or letter
m acts by R -> R o m^-1.  tau_h: (l, R) -> (l + h, c'_h^-1 o R o c_h^-1).  In character chi a term on array l in cell C
counts chi(l) on the array-0 cell sigma(c'_l) C (P_chi tau*_l = chi(l) P_chi).
A unit needs, for every character chi outside E = {eigencharacters of the shared moves}:
    e_K = b(chi) v_chi(top rep) + sum r(chi) v_chi(other free reps and inputs)          (shared moves vanish in chi)
The top rep is D(content); its content is a free combination of m.y over earlier generators y and moves m; the cell
of m.y's terms depends only on the neighbour Lam_m of L_Z (30 choices), so there are 30 menu items per generator.
Items built on a FREE orbit (or the inputs) can be taken from any member tau_j, which multiplies them by chi(j): their
coefficient is free per character.  Items built on a SHARED move have one coefficient for all characters.  So:
    for all chi not in E:  Q_chi( sum_{F,k} gamma_{F,k} v_chi(m_k F) ) = mu_chi Q_chi(e_K),   mu_chi != 0,
Q_chi = projection orthogonal to the per-character free span.  Linear in (gamma, mu); decided by SVD with a rank gap
check (entries are sums of roots of unity).
"""
import sys, os, random, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sp8 import *

class Eng:
    def __init__(self, name, P=None, seed=1):
        self.name = name
        self.P = P or Polar()
        GR = groups()
        c, cp = GR[name]
        self.c = c; self.cp = cp
        self.ci = [inv(x) for x in c]; self.cpi = [inv(x) for x in cp]
        self.W, self.add, self.neg = hindex_table(name, c)
        self.chars = characters(name)
        self.pcp = [self.P.perm(x) for x in cp]
        self.iZ = self.P.idx[LZ]
        self.iK = self.P.idx[img_mask(table(inv(sigma(G.K))), LZ)]
        self.rng = random.Random(seed)
        self.nbZ = self.P.nb[self.iZ]
        self.nbB = [basis_of(mask_elems(self.P.L[k])) for k in self.nbZ]
        pc = [self.P.perm(x) for x in c]
        self.fix = [j for j in self.nbZ if all(p[j] == j for p in pc)]
        self.Dsig = [sigma(G.D[z]) for z in range(1, 16)]
        self._letters = [sigma(G.rand_engine(self.rng)) for _ in range(400)]
        # move elements, sorted by the neighbour cell Lam_m = m^-1 L_Z they produce
        self.moves = {}
        for _ in range(3000):
            m = mul(self.rng.choice(self.Dsig), self.rng.choice(self._letters))
            k = self.P.idx[span_mask([app(inv(m), 1 << (4 + j)) for j in range(4)])]
            self.moves.setdefault(k, []).append(m)
        assert set(self.moves) == set(self.nbZ), 'not every neighbour reached'
    # ------------------------------------------------------------ terms
    def cell(self, R):
        return self.P.idx[span_mask([R[4 + j] for j in range(4)])]
    def inp(self, l):
        return {(l, IDM): 1.0 + 0j}
    def lmul(self, m, u):
        mi = inv(m)
        out = {}
        for (l, R), c in u.items():
            k = (l, mul(R, mi)); out[k] = out.get(k, 0) + c
        return {k: v for k, v in out.items() if abs(v) > 1e-12}
    def tau(self, h, u):
        out = {}
        for (l, R), c in u.items():
            k = (self.add(l, h), mul(self.cpi[h], mul(R, self.ci[h]))); out[k] = out.get(k, 0) + c
        return out
    def add_gen(self, *us, coefs=None):
        out = {}
        for i, u in enumerate(us):
            a = 1 if coefs is None else coefs[i]
            for k, v in u.items(): out[k] = out.get(k, 0) + a * v
        return {k: v for k, v in out.items() if abs(v) > 1e-12}
    def shared(self, m, y, psi=0):
        """exact (psi-twisted) orbit sum: m . sum_h psi(h)^-1 tau_h(y); absorbed iff Lam_m is H-fixed"""
        ch = self.chars[psi]
        acc = {}
        for h in range(self.W):
            th = self.tau(h, y)
            for k, v in th.items(): acc[k] = acc.get(k, 0) + np.conj(ch[h]) * v
        return self.lmul(m, {k: v for k, v in acc.items() if abs(v) > 1e-12})
    def move_to(self, k):
        return self.rng.choice(self.moves[k])
    def rand_move(self):
        return self.move_to(self.rng.choice(self.nbZ))
    def absorbed_move(self):
        return self.move_to(self.rng.choice(self.fix))
    # ------------------------------------------------------------ distances
    def dist(self, R):
        return int(self.P.dZ[self.cell(R)])
    def maxdist(self, u):
        return max([self.dist(R) for (l, R) in u] + [-1])
    def gscore(self, u):
        """geodesic progress: max over terms of d(input, cell) - d(cell, kernel) (4 at the kernel)"""
        if not u:
            return -99
        return max(int(self.P.dZ[c]) - int(self.P.D[c, self.iK]) for c in (self.cell(R) for (l, R) in u))
    def progress_move(self, u, tries=60):
        """a random move maximising the geodesic progress of u towards the kernel"""
        best = None; bs = -99
        for _ in range(tries):
            m = self.rand_move()
            sc = self.gscore(self.lmul(m, u)) + self.rng.random() * 0.5
            if sc > bs:
                bs, best = sc, m
        return best
    # ------------------------------------------------------------ character vectors
    def vchi(self, u, chi):
        v = {}
        ch = self.chars[chi]
        for (l, R), c in u.items():
            C = int(self.pcp[l][self.cell(R)])
            v[C] = v.get(C, 0) + c * ch[l]
        return v
    def menu_cells(self, u):
        """for each of the 30 neighbours Lam_m: list of (l, cell of m.term, coef)"""
        out = []
        terms = [(l, R, c) for (l, R), c in u.items()]
        for B in self.nbB:
            out.append([(l, self.P.idx[span_mask([app(R, b) for b in B])], c) for (l, R, c) in terms])
        return out
    def menu_vchi(self, mc, chi):
        ch = self.chars[chi]
        v = {}
        for (l, C, c) in mc:
            C2 = int(self.pcp[l][C]); v[C2] = v.get(C2, 0) + c * ch[l]
        return v
    # ------------------------------------------------------------ the decision
    def decide(self, free_gens, shared_gens, E, tol=1e-8):
        """free_gens: list of free-orbit reps (inputs included) visible to the top; shared_gens: list of shared moves.
        Returns (passes, info): passes iff some gamma, mu solve every chi not in E with all mu_chi != 0."""
        chis = [k for k in range(self.W) if k not in E]
        fmenus = [self.menu_cells(u) for u in free_gens]
        smenus = [self.menu_cells(F) for F in shared_gens]
        ncol = 30 * len(shared_gens)
        blocks = []; rhs = []; info = {}
        for chi in chis:
            vecs = []
            for u, mu_ in zip(free_gens, fmenus):
                vecs.append(self.vchi(u, chi))
                for mc in mu_:
                    vecs.append(self.menu_vchi(mc, chi))
            svecs = [self.menu_vchi(mc, chi) for sm in smenus for mc in sm]
            cells = sorted(set(k for v in vecs + svecs for k in v) | {self.iK})
            ci = {k: i for i, k in enumerate(cells)}
            A = np.zeros((len(cells), len(vecs)), dtype=complex)
            for j, v in enumerate(vecs):
                for k, x in v.items(): A[ci[k], j] += x
            Bm = np.zeros((len(cells), ncol), dtype=complex)
            for j, v in enumerate(svecs):
                for k, x in v.items(): Bm[ci[k], j] += x
            e = np.zeros(len(cells), dtype=complex); e[ci[self.iK]] = 1
            # orthonormal basis of span(A)
            if A.shape[1]:
                U, s, _ = np.linalg.svd(A, full_matrices=False)
                r = int((s > tol * max(1, s[0])).sum())
                U = U[:, :r]
            else:
                U = np.zeros((len(cells), 0))
            Q = lambda X: X - U @ (U.conj().T @ X)
            qe = Q(e)
            info[chi] = {'free_rank': U.shape[1], 'Qe': float(np.linalg.norm(qe))}
            blocks.append((Q(Bm), qe))
        # joint system  [QB_chi | -qe_chi e_chi] [gamma; mu] = 0
        nchi = len(chis)
        rows = sum(b.shape[0] for b, _ in blocks)
        M = np.zeros((rows, ncol + nchi), dtype=complex)
        r0 = 0
        for i, (QB, qe) in enumerate(blocks):
            M[r0:r0 + QB.shape[0], :ncol] = QB
            M[r0:r0 + QB.shape[0], ncol + i] = -qe
            r0 += QB.shape[0]
        # free chis (Qe = 0) are served by the free span alone: mu can be anything; their column is 0
        if M.shape[1] == 0:
            return True, info
        _, s, Vh = np.linalg.svd(M, full_matrices=True)
        rank = int((s > tol * max(1, s[0] if s.size else 1)).sum())
        null = Vh[rank:].conj().T            # basis of the null space
        gap = (s[rank - 1] if rank else 0, s[rank] if rank < len(s) else 0)
        mu_ok = []
        for i in range(nchi):
            if info[chis[i]]['Qe'] < tol:
                mu_ok.append(True)          # chi served without the shared moves
            else:
                mu_ok.append(null.shape[1] > 0 and np.max(np.abs(null[ncol + i, :])) > 1e-6)
        info['null_dim'] = null.shape[1]; info['sv_gap'] = gap; info['mu_ok'] = dict(zip(chis, mu_ok))
        return all(mu_ok), info
