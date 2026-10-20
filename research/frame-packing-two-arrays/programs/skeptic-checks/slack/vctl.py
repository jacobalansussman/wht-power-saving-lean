#!/usr/bin/env python3
"""v-slack own exact checks in the group algebra C[G] (contents = formal sums of group elements).

(A) o3-build's order-3 control word at block 3 (C6O3; same shape as the 8-move word at block 4):
      r1 := x1 + x2 ; r2 := x1 - s x2 ; K on both (2m moves) ; r2 := (1 - s + s^2)/2 (r1 - r2) ; r1 := r1 - r2.
    Check it is an exact formal identity, count moves, and test whether ANY gauge makes it coherent
    (row twists of the two rows after step 1 are e and s; coherence would need s in Paulis*scalars).
(B) Parity: an exact 1-array unit at block 1 with 2 moves (a lateral detour through the triangle of the dual
    polar graph), i.e. a unit with odd move deficit.  If it exists, appending it to any unit flips the parity
    of the deficit, so parity of W m - M is not a law of engine units.
"""
from fractions import Fraction as F
import itertools, json, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from vframes import matmul, ident, dirz, Xz, Zz, phase, addrmap, ONE, ZERO, I, cmul, cadd, conjT, sp_action

def key(M): return tuple(tuple(r) for r in M)

def fs_left(g, fs):
    """left-multiply formal sum fs (dict key->coef) by group element g (matrix)."""
    out = {}
    for k, c in fs.items():
        nk = key(matmul(g, [list(r) for r in k]))
        out[nk] = cadd(out.get(nk, ZERO), c)
    return {k: c for k, c in out.items() if c != ZERO}

def fs_add(a, b, cb=ONE):
    out = dict(a)
    for k, c in b.items():
        out[k] = cadd(out.get(k, ZERO), cmul(cb, c))
    return {k: c for k, c in out.items() if c != ZERO}

def fs_scale(a, c): return {k: cmul(c, v) for k, v in a.items() if cmul(c, v) != ZERO}

class Reg:
    """content of one register: list over inputs of formal sums."""
    def __init__(self, parts): self.p = parts
    def left(self, g): return Reg([fs_left(g, x) for x in self.p])
    def add(self, o, c=ONE): return Reg([fs_add(a, b, c) for a, b in zip(self.p, o.p)])
    def scale(self, c): return Reg([fs_scale(a, c) for a in self.p])

def control_word(m):
    n = 1 << m
    e = key(ident(n))
    # s = coordinate 3-cycle on the first 3 address bits: x = (x0,x1,x2) -> (x2,x0,x1) (columns of A)
    cols = [1 << c for c in range(m)]
    cols[0], cols[1], cols[2] = (1 << 1), (1 << 2), (1 << 0)
    s = addrmap(m, cols); s2 = matmul(s, s); s3 = matmul(s2, s)
    assert key(s3) == e
    Kmat = ident(n)
    for k in range(m): Kmat = matmul(Kmat, dirz(m, 1 << k))
    assert key(matmul(s, Kmat)) == key(matmul(Kmat, s)), "s must commute with K"
    x1 = Reg([{e: ONE}, {}]); x2 = Reg([{}, {e: ONE}])
    r1 = x1.add(x2)
    r2 = x1.add(x2.left(s), (F(-1), F(0)))
    rows_after_step1 = {'r1': [list(p.keys()) == [e] for p in r1.p], 'r2_x2_is_s': list(r2.p[1].keys()) == [key(s)]}
    moves = 0
    for k in range(m):
        d = dirz(m, 1 << k)
        r1 = r1.left(d); moves += 1
        r2 = r2.left(d); moves += 1
    diff = r1.add(r2, (F(-1), F(0)))
    half = (F(1, 2), F(0)); mhalf = (F(-1, 2), F(0))
    new2 = diff.scale(half).add(diff.left(s).scale(mhalf)).add(diff.left(s2).scale(half))
    r2 = new2
    r1 = r1.add(r2, (F(-1), F(0)))
    Kk = key(Kmat)
    ok1 = r1.p[0] == {Kk: ONE} and r1.p[1] == {}
    ok2 = r2.p[1] == {Kk: ONE} and r2.p[0] == {}
    # coherence: need a gauge h1, h2 with (row1) e h1 = e h2 and (row2) e h1 = s h2 modulo Paulis*scalars,
    # i.e. s in Paulis*scalars.  Test: does s act trivially on Pauli labels?
    act = sp_action(m, s)
    s_trivial_on_paulis = all(act((a, b)) == (a, b) for a in range(n) for b in range(n))
    return {'m': m, 'formal_identity_r1=Kx1': ok1, 'formal_identity_r2=Kx2': ok2, 'moves': moves,
            'Wm': 2 * m, 'deficit': 2 * m - moves, 'rows_after_step1': rows_after_step1,
            's_in_Paulis*scalars (needed for any coherent gauge)': s_trivial_on_paulis,
            'terms_in_r2_before_K': [len(p) for p in Reg([{e: ONE}, {key(s): (F(-1), F(0))}]).p]}

def odd_unit_block1():
    m = 1; n = 2
    d = dirz(m, 1)
    S = phase(m, 1)
    # the free group at block 1 modulo nothing: monomial matrices with entries i^k
    free = []
    units = [ONE, I, (F(-1), F(0)), (F(0), F(-1))]
    for perm in [(0, 1), (1, 0)]:
        for u0 in units:
            for u1 in units:
                M = [[ZERO, ZERO], [ZERO, ZERO]]
                M[perm[0]][0] = u0; M[perm[1]][1] = u1
                free.append(M)
    g = matmul(d, matmul(S, d))           # two paid letters with a free letter between
    target = d                            # K at block 1 is one paid letter
    found = []
    for b in free:
        for b2 in free:
            h = matmul(b, matmul(g, b2))  # b g b2, all letters
            # is h = c * target with c a scalar (any)?
            # ratio test
            c = None; ok = True
            for i in range(n):
                for j in range(n):
                    t = target[i][j]; x = h[i][j]
                    if t == ZERO and x != ZERO: ok = False
                    if t != ZERO:
                        # c = x / t  (exact complex division)
                        den = t[0] * t[0] + t[1] * t[1]
                        q = ((x[0] * t[0] + x[1] * t[1]) / den, (x[1] * t[0] - x[0] * t[1]) / den)
                        if c is None: c = q
                        elif c != q: ok = False
            if ok and c is not None:
                found.append({'c': [str(c[0]), str(c[1])]})
                if len(found) >= 3: break
        if found: break
    return {'2-move word b (dir S dir) b2 equals c * dir for some free b, b2': bool(found),
            'example_scalars': found[:3],
            'consequence': '1-array unit at block 1 with M = 2 = Wm + 1 (odd deficit -1); exact as a group identity'}

if __name__ == '__main__':
    res = {'control_word_block3': control_word(3), 'control_word_block4': control_word(4),
           'odd_unit_block1': odd_unit_block1()}
    for k, v in res.items():
        print(k)
        for a, b in v.items(): print('   ', a, ':', b)
    json.dump(res, open(__file__.replace('.py', '.json'), 'w'), indent=1, default=str)
