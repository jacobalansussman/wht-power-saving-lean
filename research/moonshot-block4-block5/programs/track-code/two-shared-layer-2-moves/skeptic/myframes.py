"""v-w5-twoshared: independent frame calculus in F2^8 (own code, does not import sp.py).
Vectors: ints, bits 0-3 = x, bits 4-7 = z.  omega(u,v) = x_u.z_v + z_u.x_v.
Lagrangians = 4-dim totally isotropic subspaces, stored as 256-bit membership masks.
d(L,M) = 4 - dim(L cap M).  D = transvection by e_x1: v -> v + z_1(v) e_x1 (x_1 += z_1).
Checks:
  (a) number of Lagrangians = 2295, neighbours of L0 = 30, dual polar graph distances consistent;
  (F4) every Q_U (U a 2-dim subspace of L0) is gated, and the stated corollary for N ~ L0, d(M,L0) = 2;
  (F1/F2) for every neighbour N = Lag(y) of L0 with d(N, L_D) = 2: U = L0 cap S_D(N) is 2-dim, inside H, inside N;
  (D4) no tuple (N_e, N_j, N_h, I_j) with the hypotheses of Lemma D4 survives both readings at layer <= 2.
"""
import itertools, sys

def omega(u, v):
    return (bin((u & 15) & (v >> 4)).count('1') + bin((u >> 4) & (v & 15)).count('1')) & 1

def span(vs):
    s = {0}
    for v in vs:
        if v not in s:
            s |= {w ^ v for w in s}
    return frozenset(s)

def mask(S):
    m = 0
    for v in S:
        m |= 1 << v
    return m

def dimof(m):
    c = bin(m).count('1')
    return c.bit_length() - 1

def dist(a, b):
    return 4 - dimof(a & b)

def elems(m):
    return [v for v in range(256) if (m >> v) & 1]

L0 = mask(span([1 << 4, 1 << 5, 1 << 6, 1 << 7]))

def SD(v):
    return v ^ ((v >> 4) & 1)          # x_1 += z_1

def img(f, m):
    return mask([f(v) for v in elems(m)])

def perp(vs):
    return [w for w in range(256) if all(omega(w, v) == 0 for v in vs)]

# enumerate Lagrangians by BFS over neighbours
def neighbours(m):
    E = elems(m)
    out = []
    # hyperplanes of L: kernels of nonzero functionals; pick 3-dim subspaces
    basis = []
    for v in E:
        if v and v not in span(basis):
            basis.append(v)
    hyps = set()
    for c in range(1, 16):
        # functional f(sum a_i b_i) = <c,a>; kernel
        ker = [v for v in E if True]
        # compute coordinates by brute force
        ker = []
        for a in range(16):
            v = 0
            for i in range(4):
                if (a >> i) & 1: v ^= basis[i]
            if bin(a & c).count('1') % 2 == 0:
                ker.append(v)
        hyps.add(frozenset(ker))
    for P in hyps:
        Pp = perp(list(P))
        # Lagrangians containing P: P + <w>, w in P^perp \ P (isotropic automatically, 1 vector)
        seen = set()
        for w in Pp:
            if w in P: continue
            L = mask(span(list(P) + [w]))
            if L != m and L not in seen:
                seen.add(L); out.append(L)
    return out

def all_lagrangians():
    seen = {L0}; frontier = [L0]
    adj = {}
    while frontier:
        nf = []
        for m in frontier:
            nb = neighbours(m)
            adj[m] = nb
            for n in nb:
                if n not in seen:
                    seen.add(n); nf.append(n)
        frontier = nf
    return sorted(seen), adj

if __name__ == '__main__':
    LAG, adj = all_lagrangians()
    print('Lagrangians:', len(LAG), ' neighbours of L0:', len(adj[L0]))
    assert all(len(adj[m]) == 30 for m in LAG)
    # distance consistency: graph distance (BFS) equals 4 - dim intersection, from L0 and from a few others
    def bfs(src):
        dd = {src: 0}; fr = [src]
        while fr:
            nf = []
            for m in fr:
                for n in adj[m]:
                    if n not in dd:
                        dd[n] = dd[m] + 1; nf.append(n)
            fr = nf
        return dd
    for src in LAG[::300]:
        dd = bfs(src)
        assert all(dd[m] == dist(src, m) for m in LAG)
    print('graph distance = codim of intersection: OK')
    LD = img(SD, L0)
    H = L0 & LD
    print('d(L0,L_D) =', dist(L0, LD), ' dim H =', dimof(H))
    assert all(SD(v) == v for v in elems(H))
    # S_D is symplectic and an involution
    assert all(SD(SD(v)) == v for v in range(256))
    assert all(omega(SD(u), SD(v)) == omega(u, v) for u in range(256) for v in range(0, 256, 7))
    # (F4) gatedness of Q_U for every 2-dim U in L0
    L0e = elems(L0)
    Us = set()
    for a, b in itertools.combinations([v for v in L0e if v], 2):
        Us.add(mask(span([a, b])))
    print('2-dim subspaces of L0:', len(Us))
    gated = True; cor = True
    for U in Us:
        Q = [m for m in LAG if (m & U) == U]
        assert len(Q) == 15, len(Q)   # quad: points of GQ(2,2) = 15
        for X in LAG:
            ds = [dist(X, y) for y in Q]
            mn = min(ds)
            pis = [y for y, d in zip(Q, ds) if d == mn]
            if len(pis) != 1: gated = False; break
            pi = pis[0]
            if any(dist(X, y) != mn + dist(pi, y) for y in Q): gated = False; break
        if not gated: break
        # corollary: N ~ L0, M with d(M,L0)=2 and L0 cap M = U: d(N,M) <= 2 iff U <= N, else 3
        for N in adj[L0]:
            for M in Q:
                if dist(M, L0) != 2 or (M & L0) != U: continue
                dnm = dist(N, M)
                if ((N & U) == U) != (dnm <= 2): cor = False
                if (N & U) != U and dnm != 3: cor = False
    print('(F4) every Q_U gated:', gated, '  corollary for neighbours:', cor)
    # (F1/F2): top frames N = Lag(y) are neighbours of L0 (y = D m, m an engine letter: S_m^-1 L_D, all 30 occur)
    tops = [N for N in adj[L0] if dist(N, LD) == 2]
    print('neighbours N of L0 with d(N, L_D) = 2 (F1):', len(tops))
    F2ok = True
    Uof = {}
    for N in tops:
        M = img(SD, N)
        U = L0 & M
        Uof[N] = U
        if not (dimof(U) == 2 and (U & H) == U and (U & N) == U): F2ok = False
    print('(F2) U = L0 cap S_D(N) is 2-dim, inside H, inside N, for all such N:', F2ok)
    # (D4) exhaustive at frame level
    NH = [N for N in adj[L0] if dist(N, LD) == 2]       # N_h with lambda(w_J) = d(S_D N_h, L0) = 2
    Ij_all = [I for I in LAG if dist(I, L0) == 3 and dist(I, LD) == 2]
    print('candidates I(q_j): d(I,L0)=3 and d(I,L_D)=2:', len(Ij_all))
    surv_pairs = 0; caseA = 0; caseB = 0; badB = 0
    for Ne in tops:
        Me = img(SD, Ne)
        for Nj in tops:
            Mj = img(SD, Nj)
            # junk/route hypotheses need some I with d(Nj,I)=4 and d(Ne,I)=2
            Is = [I for I in Ij_all if dist(Nj, I) == 4 and dist(Ne, I) == 2]
            if not Is: continue
            for Nh in NH:
                if dist(Me, Nh) <= 2 and dist(Mj, Nh) <= 2:
                    surv_pairs += 1
                    if Uof[Ne] == Uof[Nj]: caseB += 1
                    else: caseA += 1
            if Uof[Ne] == Uof[Nj]: badB += 1
    print('D4: (N_e, N_j) pairs admitting route+junk I:', 'computed above')
    print('D4: surviving (N_e, N_j, N_h) with both readings <= 2 and an admissible I(q_j):', surv_pairs,
          ' [case U_e != U_j:', caseA, ', case U_e = U_j:', caseB, ']')
    print('D4 case U_e = U_j: (N_e, N_j) pairs with U_e = U_j that admit an I(q_j):', badB)
    # also report: without the d(S_D N_h, L0) = 2 restriction, how many survive (control: the restriction matters)
    ctl = 0
    for Ne in tops:
        Me = img(SD, Ne)
        for Nj in tops:
            if not [I for I in Ij_all if dist(Nj, I) == 4 and dist(Ne, I) == 2]: continue
            Mj = img(SD, Nj)
            for Nh in adj[L0]:
                if dist(Me, Nh) <= 2 and dist(Mj, Nh) <= 2: ctl += 1
    print('control (N_h any neighbour, lambda(w_J) unrestricted): survivors', ctl)
