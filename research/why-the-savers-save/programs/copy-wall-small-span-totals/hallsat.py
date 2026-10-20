"""hallsat.py (d-R7): SAT model for 'frames of given dims whose read masks pass octad Hall on P112'.
Frame k = {x in F_2^9 : h_{k,j}.x = 0, j < 9 - d_k}; dual vectors h_{k,j} free 9-bit vectors, all nonzero
combinations nonzero (so dim F_k = d_k exactly).  y[k,p] <-> port p in F_k.  For every octad (1260) a perfect
matching of its 8 ports into distinct frames: m[O,p,k] -> y[k,p]; each port one frame; each frame <= 1 port.
Symmetry: frames of equal dim ordered by their first dual vector (lex, weak).  usage: hallsat.py dims seconds"""
import sys, time, itertools; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
dims = [int(x) for x in sys.argv[1].split(',')]; SECS = float(sys.argv[2])
K = len(dims); vp = IDPool(); cl = []
H = {(k, j, i): vp.id(('h', k, j, i)) for k in range(K) for j in range(9 - dims[k]) for i in range(9)}
def xor_var(lits):
    """fresh var equal to XOR of lits (Tseitin chain)"""
    cur = lits[0]
    for l in lits[1:]:
        z = vp.id(('x', cur, l))
        cl.extend([[-z, cur, l], [-z, -cur, -l], [z, -cur, l], [z, cur, -l]]); cur = z
    return cur
# independence: every nonzero combination of the dual vectors of a frame is nonzero
for k in range(K):
    m = 9 - dims[k]
    for sub in range(1, 1 << m):
        js = [j for j in range(m) if sub >> j & 1]
        bits = [H[k, js[0], i] if len(js) == 1 else xor_var([H[k, j, i] for j in js]) for i in range(9)]
        cl.append(bits)
# y[k,p] <-> AND_j (h_{k,j} . u_p == 0)
Y = {}
for k in range(K):
    for p in range(NP):
        sup = [i for i in range(9) if LAB[p] >> i & 1]
        zs = [xor_var([H[k, j, i] for i in sup]) for j in range(9 - dims[k])]   # z = h.u_p (1 = nonzero)
        y = vp.id(('y', k, p)); Y[k, p] = y
        for z in zs: cl.append([-y, -z])
        cl.append([y] + zs)
# octad matchings
for oi, oc in enumerate(OCTP):
    for p in oc:
        ms = [vp.id(('m', oi, p, k)) for k in range(K)]
        cl.append(ms)
        for k in range(K): cl.append([-ms[k], Y[k, p]])
    for k in range(K):
        lits = [vp.id(('m', oi, p, k)) for p in oc]
        enc = CardEnc.atmost(lits=lits, bound=1, vpool=vp, encoding=EncType.seqcounter)
        cl.extend(enc.clauses)
# weak symmetry breaking: equal-dim frames ordered by the top bit pattern of their first dual vector (lex on bits)
def lex_le(a, b):
    """a <= b lexicographically (bit 8 most significant); a, b lists of vars"""
    eq_prev = None
    for i in range(8, -1, -1):
        ai, bi = a[i], b[i]
        if eq_prev is None: cl.append([-ai, bi])
        else: cl.append([-eq_prev, -ai, bi])
        e = vp.id(('eq', ai, bi, eq_prev))
        # e <-> eq_prev & (ai == bi)
        base = [] if eq_prev is None else [eq_prev]
        cl.append([-e, -ai, bi]); cl.append([-e, ai, -bi])
        if eq_prev is not None: cl.append([-e, eq_prev])
        cl.append([e] + [-x for x in base] + [ai, bi]); cl.append([e] + [-x for x in base] + [-ai, -bi])
        eq_prev = e
# canonical basis of each 2-dim dual {a, b, a^b}: a <= b <= a^b (lex), i.e. a, b are its two smallest vectors
for k in range(K):
    if dims[k] == 7:
        a = [H[k, 0, i] for i in range(9)]; b = [H[k, 1, i] for i in range(9)]
        c = [xor_var([a[i], b[i]]) for i in range(9)]
        lex_le(a, b); lex_le(b, c)
for k in range(K - 1):
    if dims[k] == dims[k + 1]:
        lex_le([H[k, 0, i] for i in range(9)], [H[k + 1, 0, i] for i in range(9)])
print('dims', dims, 'vars', vp.top, 'clauses', len(cl), flush=True)
t0 = time.time()
with Solver(name=__import__('os').environ.get('SATS', 'g4'), bootstrap_with=cl) as s:
    import threading
    timer = threading.Timer(SECS, lambda: s.interrupt()); timer.start()
    r = s.solve_limited(expect_interrupt=True)
    timer.cancel()
    print('result', r, '%.1fs' % (time.time() - t0), flush=True)
    if r:
        mdl = set(l for l in s.get_model() if l > 0)
        fr = []
        for k in range(K):
            hs = [sum(1 << i for i in range(9) if H[k, j, i] in mdl) for j in range(9 - dims[k])]
            fr.append(hs)
        Sm = np.array([[all(par(h & LAB[p]) == 0 for h in hs) for p in range(NP)] for hs in fr])
        print('frames (dual vectors):', fr, 'ports', Sm.sum(1).tolist(), 'Hall deficit check', int(hall_deficit(Sm).clip(0).sum()), flush=True)
        np.save(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'out/hallsat_%s.npy') % '-'.join(map(str, dims)), np.array(fr, dtype=object), allow_pickle=True)
