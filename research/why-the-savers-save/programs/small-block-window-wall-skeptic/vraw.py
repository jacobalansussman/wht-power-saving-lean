#!/usr/bin/env python3
"""v-an2-smallnet / vraw.py (own code; numpy from our own venv bound-check/g-p112-decoder/venv).
Independent raw-array check of the frame calculus and of a two-pair shared-bridge call with side information.
Phi_U is built here as a PRODUCT of rank-one kernels Phi_b = J diag(i^{wt(b) (b.x)}) J over an F_2-orthonormal basis
(an2-smallnet builds one diagonal from the projection weight); basis independence is checked.
Network simulator: every array holds a physical complex vector and a frame (a set of F_2 vectors = subspace);
climb/descend apply Phi_new Phi_old^{-1} and are paid |dim change|; a gate needs equal frames and is free."""
import itertools, json, sys
import numpy as np
rng = np.random.default_rng(20261015)
pc = lambda x: bin(x).count('1')
dt = lambda x, y: pc(x & y) & 1

def J(m):
    n = 1 << m
    return np.array([[(-1) ** dt(x, y) for y in range(n)] for x in range(n)], dtype=complex) / np.sqrt(n)

def span(B):
    S = {0}
    for b in B: S |= {s ^ b for s in S}
    return frozenset(S)

def orthonormal_basis(S, m):
    """an F_2-orthonormal basis of subspace S (set of ints), or None."""
    S = set(S); basis = []
    rest = set(S)
    while len(span(basis)) < len(S):
        cand = [v for v in S if v and pc(v) % 2 == 1 and all(dt(v, b) == 0 for b in basis) and v not in span(basis)]
        if not cand: return None
        basis.append(min(cand))
    return basis

def Phi(S, m, Jm, basis=None):
    if basis is None: basis = orthonormal_basis(S, m)
    assert basis is not None, 'no orthonormal basis'
    n = 1 << m; M = np.eye(n, dtype=complex)
    for b in basis:
        d = np.array([1j ** ((pc(b) * dt(b, x)) % 4) for x in range(n)])
        M = M @ (Jm @ np.diag(d) @ Jm)
    return M

def Cm(m):
    C = 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]]); M = np.eye(1, dtype=complex)
    for _ in range(m): M = np.kron(M, C)
    return M

def R(u, m):
    n = 1 << m; M = np.zeros((n, n), dtype=complex)
    for x in range(n): M[x ^ u, x] = 1
    return M

def perp(S, m): return frozenset(v for v in range(1 << m) if all(dt(v, s) == 0 for s in S))

def identities(m):
    Jm = J(m); n = 1 << m; ok = {}
    full = frozenset(range(n))
    ok['Phi0=I'] = np.allclose(Phi(frozenset({0}), m, Jm), np.eye(n))
    ok['Phi_full=C^m'] = np.allclose(Phi(full, m, Jm), Cm(m))
    res = []
    for u in range(1, n):
        if pc(u) % 2 == 0: continue
        U = span([u]); P = perp(U, m)
        if orthonormal_basis(P, m) is None: continue
        A = Phi(U, m, Jm)
        res.append(np.allclose(A @ A, R(u, m)) and np.allclose(Phi(P, m, Jm) @ np.linalg.inv(A), Cm(m) @ R(u, m)))
    ok['Phi_u^2=R_u and Phi_uperp Phi_u^-1 = C^m R_u (all odd u with orthonormal complement)'] = all(res) and len(res) > 0
    # basis independence and multiplicativity on orthogonal sums, random subspaces
    bi, mu = [], []
    for _ in range(20):
        k = rng.integers(1, m + 1); vecs = list(rng.integers(1, n, size=k))
        S = span(vecs); b1 = orthonormal_basis(S, m)
        if b1 is None: continue
        # another orthonormal basis: random orthogonal transform by trying other choices
        allb = [list(p) for p in itertools.permutations([v for v in S if v and pc(v) % 2], len(b1))
                if all(dt(p[i], p[j]) == (1 if i == j else 0) for i in range(len(p)) for j in range(len(p)))
                and span(p) == S][:5]
        for b2 in allb: bi.append(np.allclose(Phi(S, m, Jm, b1), Phi(S, m, Jm, b2)))
        if len(b1) >= 2:
            U, V = span(b1[:1]), span(b1[1:])
            mu.append(np.allclose(Phi(S, m, Jm), Phi(U, m, Jm) @ Phi(V, m, Jm)))
    ok['basis independent'] = all(bi) and len(bi) > 0
    ok['Phi_{U+V} = Phi_U Phi_V, U perp V'] = all(mu) and len(mu) > 0
    return ok

class Net:
    def __init__(s, m):
        s.m = m; s.Jm = J(m); s.phys = {}; s.fr = {}; s.paid = 0; s.desc = 0; s.cache = {}
    def ph(s, S):
        if S not in s.cache: s.cache[S] = Phi(S, s.m, s.Jm)
        return s.cache[S]
    def new(s, name, S, vec): s.phys[name] = vec.copy(); s.fr[name] = S
    def move(s, name, T):
        S = s.fr[name]
        assert S <= T or T <= S, 'not nested'
        s.phys[name] = s.ph(T) @ np.linalg.solve(s.ph(S), s.phys[name])
        d = abs(len(T).bit_length() - len(S).bit_length()); s.paid += d
        if T < S: s.desc += d
        s.fr[name] = T
    def gate(s, tgt, c, src):
        assert s.fr[tgt] == s.fr[src], f'gate {tgt} {src}: frames differ'
        s.phys[tgt] = s.phys[tgt] + c * s.phys[src]

def two_pair(m):
    """ports u_a = e_0, u_b = e_1 (orthogonal labels), one shared bridge Z = x_a + x_b, side-info copies Ka, Kb."""
    n = 1 << m; ua, ub = 1, 2; N = Net(m)
    Z0 = frozenset({0}); A = span([ua]); B = span([ub]); AB = span([ua, ub])
    xa, xb, ya, yb = (rng.normal(size=n) + 1j * rng.normal(size=n) for _ in range(4))
    N.new('Xa', A, xa); N.new('Xb', B, xb); N.new('Ya', Z0, ya); N.new('Yb', Z0, yb)
    for k in ('Z', 'Ka', 'Kb'): N.new(k, Z0, np.zeros(n, dtype=complex))
    N.move('Kb', B); N.gate('Kb', 1, 'Xb')
    N.move('Ka', A); N.gate('Ka', 1, 'Xa')
    N.move('Z', A); N.gate('Z', 1, 'Xa'); N.move('Xb', AB); N.move('Z', AB); N.gate('Z', 1, 'Xb')
    N.move('Z', Z0)                                     # the bridge's descent (rank of the label pairing = 2)
    N.gate('Ya', 1, 'Z'); N.gate('Yb', 1, 'Z')
    N.move('Ya', B); N.gate('Ya', -1, 'Kb')             # junk x_b removed; Ya's frame now holds u_b
    N.move('Yb', A); N.gate('Yb', -1, 'Ka')
    N.gate('Ka', -1, 'Xa')                              # condition A: copies erased against their sources
    N.move('Kb', AB); N.gate('Kb', -1, 'Xb')
    N.move('Xa', AB); N.move('Z', AB); N.gate('Z', -1, 'Xa'); N.gate('Z', -1, 'Xb')
    ok = {}
    ok['Ya frame holds partner label u_b'] = ub in N.fr['Ya']
    full = frozenset(range(n))
    for k in ('Xa', 'Xb'): N.move(k, full)
    N.move('Ya', perp(A, m)); N.move('Yb', perp(B, m))
    Cmm = Cm(m)
    # inputs: x given physically at its line (free); outputs: Y = Phi_{u^perp}(content) where content(x) means
    # Phi_u^{-1} x_phys.  So expected Ya = C^m R_ua xa + Phi_{ua^perp} ya.
    ok['Ya = C^m R_ua x_a + Phi(ua^perp) y_a exactly'] = np.allclose(N.phys['Ya'], Cmm @ R(ua, m) @ xa + N.ph(perp(A, m)) @ ya)
    ok['Yb = C^m R_ub x_b + Phi(ub^perp) y_b exactly'] = np.allclose(N.phys['Yb'], Cmm @ R(ub, m) @ xb + N.ph(perp(B, m)) @ yb)
    ok['Xa = C^m Phi_ua^-1 x_a (unchanged data, full frame)'] = np.allclose(N.phys['Xa'], Cmm @ np.linalg.solve(N.ph(A), xa))
    ok['scratch Z, Ka, Kb zero'] = all(np.allclose(N.phys[k], 0) for k in ('Z', 'Ka', 'Kb'))
    ok['paid'] = N.paid; ok['descended levels'] = N.desc
    return ok

def lone_pair(m, u):
    n = 1 << m; N = Net(m); U = span([u]); Z0 = frozenset({0}); full = frozenset(range(n))
    x, y = (rng.normal(size=n) + 1j * rng.normal(size=n) for _ in range(2))
    N.new('X', U, x); N.new('Y', Z0, y)
    N.move('X', Z0)                     # one descent of the data array (it acts as its own bridge)
    N.gate('X', 1, 'Y'); N.gate('Y', -1, 'X'); N.gate('X', 1, 'Y')   # free swap at a common frame: X=y... check
    # after: X = x + y - (x + y) ... compute: X=x+y; Y=y-(x+y)=-x; X=x+y+(-x)=y.  So X=y, Y=-x
    N.move('Y', perp(U, m)); N.move('X', full)
    ok = dict(paid=N.paid, two_m=2 * m,
              Y_ok=np.allclose(N.phys['Y'], -(Cm(m) @ R(u, m) @ x)), X_ok=np.allclose(N.phys['X'], Cm(m) @ y))
    return ok

if __name__ == '__main__':
    out = {}
    for m in (1, 2, 3, 4):
        out[f'identities m={m}'] = {k: bool(v) for k, v in identities(m).items()}
        print(m, out[f'identities m={m}'], flush=True)
    for m in (3, 4, 5):
        u = 1 if m != 3 else 1
        if orthonormal_basis(perp(span([u]), m), m) is None: continue
        lp = lone_pair(m, u); out[f'lone pair m={m}'] = {k: (bool(v) if isinstance(v, (bool, np.bool_)) else int(v)) for k, v in lp.items()}
        print('lone pair', m, out[f'lone pair m={m}'], flush=True)
    for m in (3, 4, 5):
        tp = two_pair(m); out[f'two pair m={m}'] = {k: (bool(v) if isinstance(v, (bool, np.bool_)) else int(v)) for k, v in tp.items()}
        print('two pair', m, out[f'two pair m={m}'], flush=True)
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
