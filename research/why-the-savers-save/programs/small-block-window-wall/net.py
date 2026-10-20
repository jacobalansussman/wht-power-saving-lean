#!/usr/bin/env python3
"""an2-smallnet / net.py (own code; numpy from our venv; python3 -I).  A small exact network simulator for the frame
calculus, run on raw complex arrays: every array carries a frame (an orthonormal F_2 basis) and a physical content
(2^m x f).  step(r, basis) moves r to a nested frame, applying Phi_new Phi_old^-1, and is paid |dim change| moves
when r holds nonzero content (a zero array changes frame for free).  gate(rs, M) mixes arrays that share ONE frame
(checked as subspaces) with a scalar matrix: free.  Two networks are checked against the target exactly:
 N2  one pair, the three-shear swap y += x, x -= y, y += x, each call through a scratch centre that descends one
     level and is erased again (mixtures y + x travel on the arrays): paid 2m + 4 = W m - 2*1 + 2*(levels dropped)
     where scratch levels count once each way; outputs C^m R_u x on Y and -C^m y on X.
 N3  two pairs with orthogonal ports e1, e2 and ONE shared centre x1 + x2 in the first call: the junk x2 at Y1 and
     x1 at Y2 is removed by side information that leaves X_b while X_b's frame is still orthogonal to the other
     port (the condition of an2-smallnet NOTE sec. 2); after the call Y_a holds y_a + x_a exactly, the centre and the
     two carriers are zero, and Y_a's frame holds its junk partner's label (the target-side window lemma)."""
import numpy as np, os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frames import walsh, phi, Cm, shift, perp_basis, pc, dot

def rank2(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)

def same_space(A, B):
    return rank2(A) == rank2(B) == rank2(list(A) + list(B))

class Net:
    def __init__(s, m, f=2, seed=0):
        s.m, s.n, s.J = m, 1 << m, walsh(m)
        s.rng = np.random.default_rng(seed); s.f = f
        s.A = {}; s.paid = 0; s.cache = {}
    def P(s, basis):
        key = tuple(sorted(basis))
        if key not in s.cache: s.cache[key] = phi(list(basis), s.m, s.J)
        return s.cache[key]
    def new(s, name, basis, content=None):
        c = np.zeros((s.n, s.f), dtype=complex) if content is None else content.copy()
        s.A[name] = [list(basis), c]
    def step(s, r, basis):
        old, c = s.A[r]
        lo, hi = (old, basis) if rank2(old) <= rank2(basis) else (basis, old)
        assert rank2(list(lo) + list(hi)) == rank2(hi), 'step between frames that are not nested'
        if np.abs(c).max() > 1e-12:
            s.paid += abs(rank2(basis) - rank2(old))
            c = s.P(basis) @ np.linalg.inv(s.P(old)) @ c
        s.A[r] = [list(basis), c]
    def gate(s, rs, M):
        F = s.A[rs[0]][0]
        for r in rs: assert same_space(s.A[r][0], F), f'gate without a common frame: {r}'
        C = [s.A[r][1] for r in rs]
        for i, r in enumerate(rs):
            s.A[r][1] = sum(M[i][j] * C[j] for j in range(len(rs)))
    def zero(s, r): return np.abs(s.A[r][1]).max() < 1e-9

def N2(m):
    u = max((v for v in range(1, 1 << m) if pc(v) % 2 and len(perp_basis([v], m)) == m - 1), key=pc)
    whole = [1 << i for i in range(m)]; up = perp_basis([u], m)
    N = Net(m, seed=m)
    x = N.rng.normal(size=(N.n, 2)) + 1j * N.rng.normal(size=(N.n, 2))
    y = N.rng.normal(size=(N.n, 2)) + 1j * N.rng.normal(size=(N.n, 2))
    N.new('X', [u], x); N.new('Y', [], y)
    for call, (src, dst, sgn) in enumerate([('X', 'Y', 1), ('Y', 'X', -1), ('X', 'Y', 1)]):
        c = f'C{call}'
        hiF, loF = [u], []
        if src == 'X':      # crossing: gather high, descend, scatter low, climb back, erase against the source
            N.new(c, hiF); N.gate([src, c], [[1, 0], [1, 1]]); N.step(c, loF)
            N.gate([dst, c], [[1, sgn], [0, 1]]); N.step(c, hiF); N.gate([src, c], [[1, 0], [-1, 1]])
        else:               # free direction: gather low, climb, scatter high, descend, erase against the source
            N.new(c, loF); N.gate([src, c], [[1, 0], [1, 1]]); N.step(c, hiF)
            N.gate([dst, c], [[1, sgn], [0, 1]]); N.step(c, loF); N.gate([src, c], [[1, 0], [-1, 1]])
        assert N.zero(c)
    N.step('X', whole); N.step('Y', up)
    okY = np.allclose(N.A['Y'][1], Cm(m) @ shift(u, m) @ x)
    okX = np.allclose(N.A['X'][1], -(Cm(m) @ y))
    return dict(m=m, paid=N.paid, expected=2 * m + 4, outputs_exact=bool(okX and okY))

def N3(m):
    e1, e2 = 1, 2
    N = Net(m, seed=100 + m)
    xs = [N.rng.normal(size=(N.n, 2)) + 1j * N.rng.normal(size=(N.n, 2)) for _ in range(2)]
    ys = [N.rng.normal(size=(N.n, 2)) + 1j * N.rng.normal(size=(N.n, 2)) for _ in range(2)]
    N.new('X1', [e1], xs[0]); N.new('X2', [e2], xs[1]); N.new('Y1', [], ys[0]); N.new('Y2', [], ys[1])
    # side information first: carriers copy x_b while X_b is still on its own line (orthogonal to the other port)
    N.new('K21', [e2]); N.gate(['X2', 'K21'], [[1, 0], [1, 1]])      # K21 = x2, for target Y1
    N.new('K12', [e1]); N.gate(['X1', 'K12'], [[1, 0], [1, 1]])      # K12 = x1, for target Y2
    # the shared centre: both sources climb to <e1, e2>, one gate, descend 2 levels to 0, scatter to both targets
    F = [e1, e2]
    N.step('X1', F); N.step('X2', F); N.new('C', F)
    N.gate(['X1', 'X2', 'C'], [[1, 0, 0], [0, 1, 0], [1, 1, 1]])      # C = x1 + x2
    N.step('C', [])
    N.gate(['Y1', 'Y2', 'C'], [[1, 0, 1], [0, 1, 1], [0, 0, 1]])      # Y_a += x1 + x2 (junk: the other x)
    # junk removal: Y1 climbs to the carrier's frame <e2> (inside e1^perp), Y2 to <e1> (inside e2^perp)
    N.step('Y1', [e2]); N.gate(['Y1', 'K21'], [[1, -1], [0, 1]])
    N.step('Y2', [e1]); N.gate(['Y2', 'K12'], [[1, -1], [0, 1]])
    # erase: carriers climb to their source's frame and subtract it; the centre climbs back and subtracts both sources
    N.step('K21', F); N.gate(['X2', 'K21'], [[1, 0], [-1, 1]])
    N.step('K12', F); N.gate(['X1', 'K12'], [[1, 0], [-1, 1]])
    N.step('C', F); N.gate(['X1', 'X2', 'C'], [[1, 0, 0], [0, 1, 0], [-1, -1, 1]])
    clean = N.zero('C') and N.zero('K21') and N.zero('K12')
    # logical check: Y_a now holds Phi_{frame}(Phi_{e_a}^{-1} x_a) + Phi_frame(y_a): compare after moving to 0
    ok = True
    for a, (Y, eb) in enumerate([('Y1', e2), ('Y2', e1)]):
        content = np.linalg.inv(N.P(N.A[Y][0])) @ N.A[Y][1]
        want = ys[a] + np.linalg.inv(N.P([[e1, e2][a]])) @ xs[a]
        ok &= np.allclose(content, want)
        ok &= rank2(N.A[Y][0] + [eb]) == rank2(N.A[Y][0])          # Y_a's frame holds its junk partner's label
    return dict(m=m, paid_in_call=N.paid, scratch_levels=dict(centre=4, carriers=2), y_plus_x_exact=bool(ok),
                scratch_clean=bool(clean))

if __name__ == '__main__':
    out = {'N2': [], 'N3': []}
    for m in (2, 3, 4, 5):
        r = N2(m); out['N2'].append(r)
        print(f"N2 one-pair three-shear swap, m={m}: paid {r['paid']} (2m+4 = {r['expected']}), outputs exact {r['outputs_exact']}")
        assert r['outputs_exact'] and r['paid'] == r['expected']
    for m in (3, 4, 5):
        r = N3(m); out['N3'].append(r)
        print(f"N3 two pairs, one shared centre, side information, m={m}: Y_a = y_a + x_a exactly {r['y_plus_x_exact']},"
              f" scratch erased {r['scratch_clean']}, paid inside the call {r['paid_in_call']}")
        assert r['y_plus_x_exact'] and r['scratch_clean']
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'net.json'), 'w'), indent=1)
