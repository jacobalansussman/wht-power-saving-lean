"""sfield.py (rm-singer): F_{2^m}, a self-dual basis, the Singer cycle in self-dual coordinates, binary quadratic
forms, the linear part a_Q, and the canonical Singer data (K_Q, J_Q) of every form.  Pure Python, exact over F_2.

Conventions (match e3-rmbuild Lemma C / w8-rmexact Theorem MT):
  * Sym_m = symmetric m x m matrices over F_2; a move in direction z is [z z^T]; the target is [I].
  * Q(x) = sum_{i<=j} q_ij x_i x_j; pairing <Q,S> = sum_{i<=j} q_ij S_ij, so <Q, zz^T> = Q(z).
  * Address map h acts by S -> h S h^T. In the induced representation Ind chi_Q = C^n (channels k, P_g|k> = |gk>),
    [S] acts on channel k by (-1)^{<Q, h_k^-1 S h_k^-T>} with h_k = sigma^k.
  * Coordinates: a SELF-DUAL basis w_1..w_m of F over F_2 (Tr(w_i w_j) = delta_ij).  Then x = sum_i Tr(x w_i) w_i,
    the dot product of coordinates is Tr(xy), and sigma = (mult by alpha) is a SYMMETRIC matrix.
"""
import itertools

POLY = {2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011}


class Field:
    def __init__(s, m):
        s.m = m; s.q = 1 << m; s.n = s.q - 1; s.poly = POLY[m]
        s.exp = [1]
        for _ in range(s.n - 1):
            s.exp.append(s.mul(s.exp[-1], 2))
        assert len(set(s.exp)) == s.n, "not primitive"
        s.log = {v: k for k, v in enumerate(s.exp)}

    def mul(s, a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a >> s.m: a ^= s.poly
        return r

    def pw(s, a, e):
        r = 1
        for _ in range(e % s.n if a else e): r = s.mul(r, a)
        return r if a else (1 if e == 0 else 0)

    def a_exp(s, k):  # alpha^k
        return s.exp[k % s.n]

    def tr(s, x):
        t = 0; y = x
        for _ in range(s.m):
            t ^= y; y = s.mul(y, y)
        assert t in (0, 1)
        return t

    def selfdual_basis(s):
        m = s.m; elems = list(range(1, s.q))
        def dfs(chosen):
            if len(chosen) == m:
                return chosen
            for w in elems:
                if w in chosen: continue
                if s.tr(s.mul(w, w)) != 1: continue
                if any(s.tr(s.mul(w, c)) for c in chosen): continue
                # linear independence over F_2
                span = {0}
                for c in chosen:
                    span |= {x ^ c for x in span}
                if w in span: continue
                r = dfs(chosen + [w])
                if r: return r
            return None
        return dfs([])


class Singer:
    """Everything in self-dual coordinates."""
    def __init__(s, m, coords='sd'):
        s.F = F = Field(m); s.m = m; s.n = F.n; s.coords = coords
        s.w = F.selfdual_basis()
        assert s.w is not None
        if coords == 'sd':   # coordinates of x: c_i = Tr(x w_i)  (sigma symmetric)
            s.coord = {x: tuple(F.tr(F.mul(x, wi)) for wi in s.w) for x in range(F.q)}
        else:                # polynomial basis 1, alpha, ..., alpha^{m-1} (the placement of e3-rmbuild / w8)
            s.coord = {x: tuple((x >> i) & 1 for i in range(m)) for x in range(F.q)}
            s.w = [1 << i for i in range(m)]
        s.elem = {v: x for x, v in s.coord.items()}
        # sigma as a matrix: column i = coords of alpha * w_i
        s.sigma = [[s.coord[F.mul(2, s.w[j])][i] for j in range(m)] for i in range(m)]
        s.idx = [(i, j) for i in range(m) for j in range(i, m)]
        s.N = len(s.idx)

    # --- matrices over F_2 as lists of lists
    @staticmethod
    def mm(A, B):
        n, k, p = len(A), len(B), len(B[0])
        return [[sum(A[i][t] & B[t][j] for t in range(k)) & 1 for j in range(p)] for i in range(n)]

    @staticmethod
    def tp(A):
        return [list(r) for r in zip(*A)]

    def mpow(s, A, e):
        m = s.m; R = [[int(i == j) for j in range(m)] for i in range(m)]
        e %= s.n
        for _ in range(e): R = s.mm(R, A)
        return R

    def pair(s, Q, S):  # Q: tuple of q_ij in s.idx order; S: symmetric matrix
        return sum(Q[b] & S[i][j] for b, (i, j) in enumerate(s.idx)) & 1

    def Qval(s, Q, x):  # x: field element
        c = s.coord[x]
        return sum(Q[b] & c[i] & c[j] for b, (i, j) in enumerate(s.idx)) & 1

    def linpart(s, Q):
        """a_Q = coefficient of x in the univariate polynomial of Q: a = sum_{y != 0} Q(y) y^{-1}."""
        F = s.F; a = 0
        for y in range(1, F.q):
            if s.Qval(Q, y): a ^= F.a_exp(-F.log[y])
        return a

    def forms(s):
        return itertools.product((0, 1), repeat=s.N)

    def KJ(s, Q, z0=1):
        """K = {k: Q(alpha^-k z0) = 1}; J = {k: <Q, sigma^-k sigma^-kT> = 1} (the target's channels)."""
        F = s.F; n = s.n
        K = frozenset(k for k in range(n) if s.Qval(Q, F.mul(F.a_exp(-k), z0)))
        J = frozenset(k for k in range(n) if s.pair(Q, s.Tk(k)))
        return K, J

    def Tk(s, k):
        if not hasattr(s, '_Tk'):
            s._Tk = {}
        k %= s.n
        if k not in s._Tk:
            hinv = s.mpow(s.sigma, -k)
            s._Tk[k] = s.mm(hinv, s.tp(hinv))
        return s._Tk[k]
