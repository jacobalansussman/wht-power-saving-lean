"""lift.py (pk2-complete): frame-game programs, their lift to explicit ENGINE WORDS, and a plain evaluator.

A frame-game program (module form, fp-build PROOF Definition 2.1): states built by
    ('in', l)                      the input l
    ('mix', [(x, p, r), ...])      sum x * p * state_r, x rational, p in P   (free)
    ('move', z, r)                 T_z state_r                               (paid)
and, for each input l, a readout [(x, p, r), ...] with sum x p state_r = Kbar e_l, found by fpa.end_solve (exact, Q).

The lift (complete.py): every free p becomes Lambda(p) = [p^] D(p^); a move becomes the free step
D(dir z) . state_r into a scratch register followed by the paid letter [dir z]; the readout becomes the free step
sum x Lambda(p) state_r followed by the free correction Q.  The resulting ENGINE WORD is stored as plain data:
free steps are lists of (coef, head matrix, u, j, src) meaning coef * [head . i^j X^a_u Z^b_u] applied to register
src; paid steps are (z, src) meaning [dir z].  The evaluator below knows nothing about Lambda: it applies group
elements g as g^(x)f, exactly as the engine does at f digit columns.  It also checks that every head lies in B
(symplectic image in P) and counts paid letters.  (Free-step terms are stored as (float coef, rational x, complex c,
head, u, j, src) with coef = x * c, so that the mod-p evaluator can use x and c exactly.)

Usage (block m, W inputs):
    from lift import *
    E = Eng(m); P = FrameProg(E.S, W)
    y = P.inp(0); r = P.mix([(1, p, y), (Fr(1, 2), q, y)]); r = P.move(z, r); ...   # p, q in P (fpa matrices)
    assert P.solve() == 'yes'                  # exact readouts by fpa.end_solve, re-verified
    word = lift(E, P)                          # explicit engine word (dict: steps, heads, outs, W)
    run_float(E, word, f, seed); run_modp(E, word, f, seed, 167772161)
"""
import numpy as np
from fractions import Fraction as Fr
from complete import *

class FrameProg:
    def __init__(s, S, W):
        s.S = S; s.W = W; s.ops = []; s.states = []
    def _add(s, op, c): s.ops.append(op); s.states.append(c); return len(s.states) - 1
    def inp(s, l): return s._add(('in', l), inp(s.S, l))
    def mix(s, terms):
        c = {}
        for x, p, r in terms: c = cadd((1, c), (Fr(x), lmul(s.S, p, s.states[r])))
        return s._add(('mix', [(Fr(x), p, r) for x, p, r in terms]), c)
    def move(s, z, r): return s._add(('move', z, r), lmul(s.S, s.S.T(z), s.states[r]))
    def moves(s): return sum(1 for o in s.ops if o[0] == 'move')
    def solve(s):
        s.read = []
        for l in range(s.W):
            st, info = end_solve(s.S, s.states, {(l, s.S.K): Fr(1)})
            if st != 'yes': return st
            sol = info['solution']
            ok, _ = verify_combo(s.S, s.states, sol, {(l, s.S.K): Fr(1)})
            assert ok
            s.read.append(sol)
        return 'yes'

def lift(E, prog, mode='full'):
    """mode 'full': Lambda = sum_k Lambda_k e_k (this track); 'f1': fp-build's f = 1 lift, D = Pexp(g^-1) unprojected"""
    if mode == 'full':
        Dof = E.D; Q = E.Qout()
    else:
        Dof = lambda M: {(u, 0): c for u, c in enumerate(E.pexp(np.linalg.inv(M))) if abs(c) > 1e-13}
        Q = {(u, 0): c for u, c in enumerate(E.pexp(E.K)) if abs(c) > 1e-13}
    I = np.eye(E.n, dtype=complex)
    heads = [I]; hid = {}
    def head(p):
        if p not in hid: hid[p] = len(heads); heads.append(E.lift_P(p))
        return hid[p]
    Dcache = {}
    def Dp(p):
        if p not in Dcache: Dcache[p] = Dof(heads[head(p)])
        return Dcache[p]
    Dz = {}
    steps = []; regmap = {}
    for i, op in enumerate(prog.ops):
        if op[0] == 'in':
            steps.append(('in', op[1]))
        elif op[0] == 'mix':
            terms = []
            for x, p, r in op[1]:
                for (u, j), c in Dp(p).items(): terms.append((complex(c) * float(x), x, c, head(p), u, j, regmap[r]))
            steps.append(('free', terms))
        else:
            z, r = op[1], op[2]
            if z not in Dz: Dz[z] = Dof(E.dirm(z))
            steps.append(('free', [(c, Fr(1), c, 0, u, j, regmap[r]) for (u, j), c in Dz[z].items()]))
            steps.append(('paid', z, len(steps) - 1))
        regmap[i] = len(steps) - 1
    outs = []
    for l in range(prog.W):
        terms = []
        for r, p, x in prog.read[l]:
            for (u, j), c in Dp(p).items(): terms.append((complex(c) * float(x), x, c, head(p), u, j, regmap[r]))
        steps.append(('free', terms)); a = len(steps) - 1
        steps.append(('free', [(c, Fr(1), c, 0, u, j, a) for (u, j), c in Q.items()]))
        outs.append(len(steps) - 1)
    return dict(steps=steps, heads=heads, outs=outs, W=prog.W)

def word_stats(E, word):
    paid = sum(1 for s in word['steps'] if s[0] == 'paid')
    free_terms = sum(len(s[1]) for s in word['steps'] if s[0] == 'free')
    for h in word['heads']: assert E.S.in_P(E.symp(h)), "a free head outside B"
    return paid, free_terms

def last_use(word):
    """index of the last step that reads each register (outputs are kept)"""
    lu = {}
    for i, st in enumerate(word['steps']):
        if st[0] == 'paid': lu[st[2]] = i
        elif st[0] == 'free':
            for t in st[1]: lu[t[-1]] = i
    for o in word['outs']: lu[o] = len(word['steps'])
    return lu

def run_float(E, word, f, seed):
    R = Raw(E, f); rn = np.random.default_rng(seed)
    xs = [rn.normal(size=R.N) + 1j * rn.normal(size=R.N) for _ in range(word['W'])]
    regs = []; lu = last_use(word)
    for i, st in enumerate(word['steps']):
        for k in [k for k, v in lu.items() if v == i - 1 and k < len(regs)]: regs[k] = None
        if st[0] == 'in': regs.append(xs[st[1]])
        elif st[0] == 'paid': regs.append(R.mat(E.dirm(st[1]), regs[st[2]]))
        else:
            groups = {}
            for (cf, x, c, h, u, j, src) in st[1]:
                groups.setdefault((h, src), {}); g = groups[(h, src)]
                g[(u, j)] = g.get((u, j), 0) + cf
            out = np.zeros(R.N, complex)
            for (h, src), Dd in groups.items():
                y = R.nfree(Dd, regs[src])
                out += y if h == 0 else R.mat(word['heads'][h], y)
            regs.append(out)
    errs = []
    for l, o in enumerate(word['outs']):
        t = R.mat(E.K, xs[l]); errs.append(float(np.abs(regs[o] - t).max() / np.abs(t).max()))
    return max(errs)

# ---------------------------------------------------------------------- exact evaluation modulo a prime p = 1 mod 4
class ModP:
    def __init__(s, p):
        s.p = p; s.iota = pow(3, (p - 1) // 4, p); assert s.iota * s.iota % p == p - 1
        s.inv2 = pow(2, p - 2, p)
    def cx(s, z):
        """exact dyadic Gaussian z -> F_p"""
        for e in range(0, 60):
            w = z * (1 << e); a, b = round(w.real), round(w.imag)
            if abs(w.real - a) < 1e-7 and abs(w.imag - b) < 1e-7:
                return (a + b * s.iota) * pow(s.inv2, e, s.p) % s.p
        raise ValueError("not dyadic: %r" % z)
    def fr(s, x): return x.numerator % s.p * pow(x.denominator % s.p, s.p - 2, s.p) % s.p
    def mat(s, M): return np.array([[s.cx(complex(v)) for v in row] for row in M], dtype=np.int64)

def run_modp(E, word, f, seed, p):
    F = ModP(p); R = Raw(E, f); rn = np.random.default_rng(seed); n = E.n
    xs = [rn.integers(0, p, size=R.N, dtype=np.int64) for _ in range(word['W'])]
    def matp(Mp, v):
        t = v.reshape((n,) * f)
        for ax in range(f):
            t = np.tensordot(Mp, t, axes=([1], [ax])) % p; t = np.moveaxis(t, 0, ax)
        return t.reshape(-1)
    def paulip(u, v):
        src, sg = R.pauli_idx(u)
        w = v[src]; return np.where(sg > 0, w, (p - w) % p)
    hp = [F.mat(h) for h in word['heads']]
    dirp = {z: F.mat(E.dirm(z)) for z in range(1, n)}
    ipow = [1, F.iota, p - 1, (p - F.iota) % p]
    regs = []; lu = last_use(word)
    for i, st in enumerate(word['steps']):
        for k in [k for k, v in lu.items() if v == i - 1 and k < len(regs)]: regs[k] = None
        if st[0] == 'in': regs.append(xs[st[1]])
        elif st[0] == 'paid': regs.append(matp(dirp[st[1]], regs[st[2]]))
        else:
            groups = {}
            for (cf, x, c, h, u, j, src) in st[1]:
                g = groups.setdefault((h, src), {})
                cc = F.fr(x) * F.cx(complex(c)) % p * ipow[(j * f) % 4] % p
                g[u] = (g.get(u, 0) + cc) % p
            out = np.zeros(R.N, dtype=np.int64)
            for (h, src), Dd in groups.items():
                y = np.zeros(R.N, dtype=np.int64)
                for u, cc in Dd.items():
                    if cc: y = (y + cc * paulip(u, regs[src])) % p
                out = (out + (y if h == 0 else matp(hp[h], y))) % p
            regs.append(out)
    Kp = F.mat(E.K)
    return all(np.array_equal(regs[o], matp(Kp, xs[l])) for l, o in enumerate(word['outs']))
