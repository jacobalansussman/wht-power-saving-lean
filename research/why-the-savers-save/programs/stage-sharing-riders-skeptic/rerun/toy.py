"""r-riders / toy.py -- the exact two-stage toy of s-wide / v-s-wide, written as formal sums.

Sub-game: block size 3, letters = rational mixing (free), shift (free), single directions (paid, 1 move each).
No address maps, no phases (the same sub-game as s-wide's handover_toy.py and v-s-wide's rider_check.py).
Fourier side: dir(z) = diag(i^chi_z(y)), shift(w) = diag((-1)^chi_w(y)), so the matrix of any word on one array is
diag(i^q) for a FRAME q: F2^3 -> Z4 (exact; no quotient by shifts anywhere).  A register is a formal sum
{(input, q): rational coefficient}; terms merge only when input AND exact q agree.  By w-engine line 8 /
w-size 2.2 this is the same as "right for every number f of digit columns".

Ports (s-wide): labels a = e1, e2, e3 (as ints 1, 2, 4).  Inputs x_T (P-bank, live array X_T) and y_T (0-bank,
live array Y_T).  Per-port route (the toy's rho): input l has LOCAL start frame s_l (x_T: chi_a, y_T: 0); a term
(l, physical q) sits at local frame q + s_l.  The hub's local rule compares LOCAL frames.
"""
import itertools
Y = range(8)
def chi(z): return tuple(bin(y & z).count('1') & 1 for y in Y)
def add(q, r): return tuple((a + b) & 3 for a, b in zip(q, r))
def sub(q, r): return tuple((a - b) & 3 for a, b in zip(q, r))
def smul(k, q): return tuple((k * a) & 3 for a in q)
Z0 = tuple(0 for _ in Y)
DIRS = list(range(1, 8))
E = [1, 2, 4]
F = add(add(chi(1), chi(2)), chi(4))                     # kernel frame (popcount)
PORTS = {a: [d for d in E if d != a] for a in E}         # D_T in climb order (s-wide)
PSI = {a: add(chi(PORTS[a][0]), chi(PORTS[a][1])) for a in E}
INPUTS = [('x', 1), ('x', 2), ('x', 4), ('y', 1), ('y', 2), ('y', 4)]
START = {('x', a): chi(a) for a in E}
START.update({('y', a): Z0 for a in E})
def phys(inp, local): return sub(local, START[inp])     # physical frame of an input sitting at a local frame
def local_of(inp, q): return add(q, START[inp])

def group():
    """all frames reachable by directions and shifts (exact Z4 functions)."""
    G = {Z0}; frontier = [Z0]
    gens = [chi(z) for z in DIRS]
    while frontier:
        nf = []
        for q in frontier:
            for g in gens:
                r = add(q, g)
                if r not in G: G.add(r); nf.append(r)
        frontier = nf
    return G

def targets(version):
    """required end contents {register: {(input, physical q): coeff}} for the two-stage toy.
    A: the full two-stage end state (stage i fwd: Y_T += X_T; stage i+1 rev: X_S -= Y_S), local frames as s-wide:
       X_S = -y_S @ F,   Y_S = (y_S + x_S) @ psi_S.
    B: s-wide's question made complete: the stage-i end state plus the coded total delivered in register Z:
       X_T = x_T @ F, Y_T = (y_T + x_T) @ psi_T, Z = sum_T (y_T + x_T) @ F.
    (all frames LOCAL; converted to physical here)."""
    T = {}
    for a in E:
        x, y = ('x', a), ('y', a)
        if version == 'A':
            T['X%d' % a] = {(y, phys(y, F)): -1}
            T['Y%d' % a] = {(y, phys(y, PSI[a])): 1, (x, phys(x, PSI[a])): 1}
        else:
            T['X%d' % a] = {(x, phys(x, F)): 1}
            T['Y%d' % a] = {(y, phys(y, PSI[a])): 1, (x, phys(x, PSI[a])): 1}
    if version == 'B':
        T['Z'] = {}
        for a in E:
            T['Z'][(('y', a), phys(('y', a), F))] = 1
            T['Z'][(('x', a), phys(('x', a), F))] = 1
    return T
