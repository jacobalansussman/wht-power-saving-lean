"""v-rm-theory: own check of rm-theory's canonical-placement claims (section 2, 3, 8.2-8.4, 8.8) at m = 3..6.
Own GF(2^m) code (primitive polys chosen here), own self-dual basis search, own form/target pairing.
Checks:
 (F1) sum f_j = 1 and sum_j f_j^(1+2^i) = 0 for i=1..m-1 (self-dual basis)
 (F2) tr Q = sum_j Q(f_j) equals Tr(b) for Q = Tr(bx) + sum_i Tr(c_i x^(1+2^i))  (random / exhaustive)
 (F4) multiplication matrices are symmetric in self-dual coordinates; sigma^k sigma^kT = mult by alpha^(2k)
 (F3/rho) Singer address group: n_odd(Q) = 2^(m-1) [b != 0]; rho = max n_odd/wt = 2 (canonical placement)
 (8.8) every nonzero b is the linear part of some (m-2)-flat form (weight 2^(m-2)), m = 4, 5
 (placement b) m = 4: for every target T congruent to I (420): rho(T) <= 3; which T give rho = 2
 (8.4) the Singer-invariant complements of Alt in Sym: count them (must be 1) at m = 3, 4, 5
"""
import itertools, random, sys
import numpy as np

PRIM = {3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011}

def run(m):
    P = PRIM[m]; q = 1 << m
    def mul(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a & q: a ^= P
        return r
    # alpha = 2 (x); check primitive
    pw = [1]
    for _ in range(q - 2): pw.append(mul(pw[-1], 2))
    assert len(set(pw)) == q - 1
    def power(a, e):
        r = 1
        for _ in range(e): r = mul(r, a)
        return r
    def tr(a):
        s = 0; x = a
        for _ in range(m): s ^= x; x = mul(x, x)
        assert s in (0, 1)
        return s
    TR = [tr(a) for a in range(q)]
    # self-dual basis: search
    elems = list(range(1, q))
    sd = None
    for cand in itertools.combinations(elems, m):
        ok = all(TR[mul(a, b)] == (1 if a == b else 0) for a in cand for b in cand)
        if ok: sd = cand; break
    assert sd is not None
    f = sd
    out = [f"m={m}: self-dual basis {f}"]
    # F1
    s = 0
    for x in f: s ^= x
    f1 = (s == 1)
    f1b = all((lambda: (lambda acc: acc)(eval_sum(i)))() == 0 for i in range(1, m)) if False else None
    ok_f1b = True
    for i in range(1, m):
        acc = 0
        for x in f: acc ^= mul(x, power(x, 1 << i))
        ok_f1b &= (acc == 0)
    out.append(f"  F1: sum f_j = 1: {f1}; sum f_j^(1+2^i) = 0 all i: {ok_f1b}")
    # coordinates
    def coord(a): return tuple(TR[mul(a, x)] for x in f)
    # element from coordinates: sum c_j f_j (since dual basis = f)
    def elem(v):
        a = 0
        for j, c in enumerate(v):
            if c: a ^= f[j]
        return a
    vecs = [tuple((x >> j) & 1 for j in range(m)) for x in range(q)]  # coordinate vectors, index = bitmask
    vidx = {v: i for i, v in enumerate(vecs)}
    elem_of = [elem(v) for v in vecs]   # element for each coordinate bitmask
    coord_idx = [vidx[coord(a)] for a in range(q)]  # bitmask of element
    assert all(coord_idx[elem_of[i]] == i for i in range(q))
    # F2: forms as truth tables on bitmask points
    gold = list(range(1, m // 2 + 1))
    def form_table(b, cs):
        tab = np.zeros(q, dtype=np.int64)
        for i in range(q):
            x = elem_of[i]
            v = TR[mul(b, x)]
            for gi, c in zip(gold, cs):
                v ^= TR[mul(c, mul(x, power(x, 1 << gi)))]
            tab[i] = v
        return tab
    rng = random.Random(1)
    nchk = 0; okF2 = True
    trials = 400 if m >= 5 else 200
    for _ in range(trials):
        b = rng.randrange(q); cs = [rng.randrange(q) for _ in gold]
        tab = form_table(b, cs)
        trq = 0
        for j in range(m): trq ^= int(tab[vidx[tuple(1 if k == j else 0 for k in range(m))]])
        okF2 &= (trq == TR[b]); nchk += 1
    out.append(f"  F2: tr Q = Tr(b) on {nchk} random forms: {okF2}")
    # matrices over F2 in coordinates: multiplication by beta: column j = coord(beta * f_j)
    def mulmat(beta):
        M = np.zeros((m, m), dtype=np.int64)
        for j in range(m):
            c = coord(mul(beta, f[j]))
            M[:, j] = c
        return M
    sym_all = all((mulmat(b) == mulmat(b).T).all() for b in range(1, q))
    sig = mulmat(2)
    okF4 = True
    Sk = np.eye(m, dtype=np.int64)
    for k in range(q - 1):
        T = (Sk @ Sk.T) % 2
        okF4 &= (T == mulmat(power(2, 2 * k))).all()
        Sk = (sig @ Sk) % 2
    out.append(f"  F4: all multiplication matrices symmetric: {sym_all}; sigma^k sigma^kT = mult(alpha^2k): {okF4}")
    # forms via symmetric pairing: Q given by upper-tri coefficient vector; Q(x) = sum_{i<=j} q_ij x_i x_j
    pairs = [(i, j) for i in range(m) for j in range(i, m)]
    nf = 1 << len(pairs)
    if m <= 5:
        # all forms tables
        X = np.array(vecs, dtype=np.int64)  # q x m
        mon = np.array([X[:, i] * X[:, j] for (i, j) in pairs]).T % 2  # q x npairs
        Qs = np.array([[(qq >> t) & 1 for t in range(len(pairs))] for qq in range(nf)], dtype=np.int64)
        tabs = (Qs @ mon.T) % 2  # nf x q
        wt = tabs[:, 1:].sum(1)
        def pairing(T):  # <Q,T> = sum q_ii T_ii + sum_{i<j} q_ij T_ij
            tv = np.array([T[i, j] for (i, j) in pairs], dtype=np.int64)
            return (Qs @ tv) % 2
        # Singer copies of target I: T_k = sigma^k sigma^kT
        Tk = []
        Sk = np.eye(m, dtype=np.int64)
        for k in range(q - 1):
            Tk.append((Sk @ Sk.T) % 2); Sk = (sig @ Sk) % 2
        nodd = sum(pairing(T) for T in Tk)
        nz = np.arange(nf) != 0
        rho = max(nodd[nz] / wt[nz])
        vals = sorted(set(nodd.tolist()))
        out.append(f"  canonical placement: n_odd values {vals}; rho = {rho}")
        # linear part b of each form: need b with Q(x) - pure part... compute: tr-parity check per copy equals Tr(b alpha^-k)
        # flats of weight 2^(m-2): their parity vectors over copies; check all nonzero 'b' patterns occur
        if m >= 4:
            flats = np.where(wt == (1 << (m - 2)))[0]
            pats = set(tuple(pairing(T)[flats_i] for T in Tk) for flats_i in flats[:0])
            # linear part via: p_k(Q) = Tr(b alpha^-k)? identify b from parity vector over k
            linpat = {}
            for bb in range(q):
                linpat[tuple(TR[mul(bb, power(2, k % (q - 1)))] for k in range(q - 1))] = bb  # convention alpha^(+k): copy k target mult(alpha^2k)
            PP = np.array([pairing(T) for T in Tk]).T  # nf x (q-1)
            okpat = True; lps = set()
            for Fi in range(nf):
                key = tuple(PP[Fi].tolist())
                if key not in linpat: okpat = False; break
                if Fi in set(flats.tolist()): lps.add(linpat[key])
            out.append(f"  F3: every form's copy-parity vector is (Tr(b alpha^-k))_k for some b: {okpat}; "
                       f"linear parts on weight-2^(m-2) forms: {len(lps - {0})} nonzero of {q-1}")
        if m == 4:
            # every target congruent to I: invertible symmetric with nonzero diagonal
            targets = []
            for bits in range(1 << len(pairs)):
                T = np.zeros((m, m), dtype=np.int64)
                for t, (i, j) in enumerate(pairs):
                    if (bits >> t) & 1: T[i, j] = T[j, i] = 1
                if T.diagonal().any() and round(abs(np.linalg.det(T))) % 2 == 1:
                    targets.append(T)
            rhos = {}
            canon = 0
            Ufield = [mulmat(b) for b in range(1, q)]
            for T in targets:
                Sk = np.eye(m, dtype=np.int64); nod = np.zeros(nf, dtype=np.int64)
                for k in range(q - 1):
                    nod += pairing((Sk @ T @ Sk.T) % 2); Sk = (sig @ Sk) % 2
                r = max(nod[nz] / wt[nz])
                rhos[r] = rhos.get(r, 0) + 1
                if abs(r - 2) < 1e-9:
                    canon += any((T == U).all() for U in Ufield)
            out.append(f"  m=4 all {len(targets)} targets: rho histogram {rhos}; rho=2 targets lying in U: {canon}")
    # 8.4 uniqueness: Singer-invariant complements of Alt in Sym, count via linear algebra over F2 (m<=5)
    if m <= 5:
        # complement = graph of a Singer-equivariant map Sym/Alt -> Alt. Sym/Alt ~ diagonal (dim m), Alt dim m(m-1)/2.
        # Equivariant maps phi: for basis d of diag, phi(sigma.d) = sigma.phi(d) mod structure. Solve linear system.
        nA = m * (m - 1) // 2
        alt_pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
        def act(S): return (sig @ S @ sig.T) % 2
        def diag_part(S): return np.diag(np.diag(S))
        def alt_vec(S):
            return np.array([S[i, j] for (i, j) in alt_pairs], dtype=np.int64)
        # unknown phi: diag basis E_jj -> alt vector (m * nA unknowns). Condition: for each j,
        # act(E_jj + phi(E_jj)) must equal D' + phi(D') where D' = diag part of act(E_jj + phi(E_jj)).
        # Linear in phi over F2: alt(act(E_jj)) + alt(act(phi_j)) = sum_k D'_kk phi_k, with D' = diag(act(E_jj)) + diag(act(phi_j)); act(alt) has zero diag.
        nun = m * nA
        rows = []; rhs = []
        def phiS(j, t):  # alt matrix for unknown index (j,t)
            S = np.zeros((m, m), dtype=np.int64); i1, i2 = alt_pairs[t]; S[i1, i2] = S[i2, i1] = 1; return S
        for j in range(m):
            E = np.zeros((m, m), dtype=np.int64); E[j, j] = 1
            AE = act(E); dE = np.diag(AE)
            # equation per alt coordinate s: alt(AE)_s + sum_t phi_{j,t} alt(act(phiS(j,t)))_s + sum_k dE_k phi_{k,s} = 0
            for s_ in range(nA):
                row = np.zeros(nun, dtype=np.int64)
                for t in range(nA):
                    row[j * nA + t] ^= alt_vec(act(phiS(j, t)))[s_]
                for k in range(m):
                    if dE[k]: row[k * nA + s_] ^= 1
                rows.append(row); rhs.append(alt_vec(AE)[s_])
        A_ = np.array(rows) % 2; b_ = np.array(rhs) % 2
        # rank over F2 of A_ and [A_|b_]
        def rk2(Mx):
            Mx = Mx.copy() % 2; r = 0
            for c in range(Mx.shape[1]):
                pv = [i for i in range(r, Mx.shape[0]) if Mx[i, c]]
                if not pv: continue
                Mx[[r, pv[0]]] = Mx[[pv[0], r]]
                for i in range(Mx.shape[0]):
                    if i != r and Mx[i, c]: Mx[i] ^= Mx[r]
                r += 1
            return r
        ra = rk2(A_); rab = rk2(np.hstack([A_, b_[:, None]]))
        out.append(f"  8.4: equivariant complements: system rank {ra} of {nun} unknowns, consistent {ra == rab}; "
                   f"number of complements = 2^{nun - ra if ra == rab else 'none'}")
    return "\n".join(out)

def eval_sum(i): return 0
for m in [int(a) for a in sys.argv[1:]] or [3, 4, 5]:
    print(run(m), flush=True)
