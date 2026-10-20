#!/usr/bin/env python3
"""cx-mirror / screen.py (own code, stdlib).  Can U8's own gate sequence carry family-specific second lines?

Exact frame-legality propagation at the FIXED gate order (pr717's rule, legal iff G <= U, as an-synth 4.1a states
it; an-build used the one-bit version for shared s):
  G(gate) = join of the start frames that reach the gate along time-respecting gate journeys (forward pass);
  U(gate) = meet of the end-frame ceilings reachable from the gate (backward pass);
  a gate is legal iff G(gate) <= U(gate); otherwise some register on it must drop (paid) before or after.
Word: U8X phase A, the scatter, phase B, as in our gcert/1 lowering (the Lean route's two dirt gates g1 and g8 are left
out: they join all y's at time 0 and all x's at the end, which two-sided frames forbid outright).  The scatter is
modelled as 8 fresh copies, one per total, that start at frame 0 after the paid drop and are read by their 56 y's in
table order (one copy per (total, reader) is also tried: the most generous scatter).  Slots start at frame 0.
Frames live in F_2^10 (window bits 0..8, s = bit 9).
Calibration: the one-sided contract (x_b from <u_b> inside the window, y_t from 0 inside u_t^perp & window) must give
0 illegal gates.  Cases: the natural mirror v_b = s + j + u_b; a symmetric twisted mirror (twist.json, best A)."""
import gzip, json, os, sys
P = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz'))
d = json.load(gzip.open(P)); v = d['v']; R = d['R']; U = d['ports']
S = 1 << 9; J9 = (1 << 9) - 1; FULL10 = (1 << 10) - 1
dot = lambda a, b: bin(a & b).count('1') & 1
def red(vs):
    basis = []
    for x in vs:
        for b in basis: x = min(x, x ^ b)
        if x:
            basis = [min(b, b ^ x) for b in basis] + [x]
    return tuple(sorted(basis, reverse=True))
def join(*sps): return red([x for sp in sps for x in sp])
def perp(sp, nb=10):
    return red([w for w in range(1, 1 << nb) if all(dot(w, x) == 0 for x in sp)]) if len(sp) < 6 else perp_fast(sp, nb)
def perp_fast(sp, nb=10):
    # nullspace of the rows of sp
    rows = list(sp); basis = []
    # solve x . r = 0 for all r: Gaussian elimination on columns
    M = [[(r >> c) & 1 for c in range(nb)] for r in rows]
    pc = []; rr = 0
    for c in range(nb):
        p = next((i for i in range(rr, len(M)) if M[i][c]), None)
        if p is None: continue
        M[rr], M[p] = M[p], M[rr]
        for i in range(len(M)):
            if i != rr and M[i][c]: M[i] = [a ^ b for a, b in zip(M[i], M[rr])]
        pc.append(c); rr += 1
    for f in range(nb):
        if f in pc: continue
        x = 1 << f
        for i, c in enumerate(pc):
            if M[i][f]: x |= 1 << c
        basis.append(x)
    return red(basis)
def meet(*sps):
    return perp_fast(join(*[perp_fast(sp) for sp in sps]))
def inside(a, b):
    bb = red(b)
    return all(red(list(bb) + [x]) == bb for x in a)
WIN = red([1 << i for i in range(9)])

def regs(g): return [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == 'out' else [])
scat = d['scat']['table']; readers = {k: [] for k in range(8)}
for t, row in enumerate(scat):
    for k, _, _ in row: readers[k].append(t)

def sequence(per_reader):
    seq = [('A', regs(g)) for g in d['A']]
    base = 2 * v + R; nc = 0
    for k in range(8):
        for t in readers[k]:
            c = base + (nc if per_reader else k)
            seq.append(('S', [c, v + t]))
            if per_reader: nc += 1
    ncopies = nc if per_reader else 8
    seq += [('B', regs(g)) for g in d['B']]
    return seq, base, ncopies

def run(start, ceil, per_reader, label):
    seq, base, nc = sequence(per_reader)
    n = base + nc
    G = [start(r) for r in range(n)]
    Gs = []
    for ph, rs in seq:
        jn = join(*[G[r] for r in rs])
        for r in rs: G[r] = jn
        Gs.append(jn)
    Uc = [ceil(r) for r in range(n)]
    Us = [None] * len(seq)
    for i in range(len(seq) - 1, -1, -1):
        ph, rs = seq[i]
        mt = meet(*[Uc[r] for r in rs])
        for r in rs: Uc[r] = mt
        Us[i] = mt
    bad = [i for i in range(len(seq)) if not inside(Gs[i], Us[i])]
    byph = {}
    for i in bad: byph[seq[i][0]] = byph.get(seq[i][0], 0) + 1
    touched = set(r for i in bad for r in seq[i][1])
    excess = sum(len(Gs[i]) - len(meet(Gs[i], Us[i])) for i in bad)
    res = dict(case=label, gates=len(seq), illegal=len(bad), by_phase=byph,
               x_on_illegal=sum(1 for r in touched if r < v), y_on_illegal=sum(1 for r in touched if v <= r < 2 * v),
               scratch_on_illegal=sum(1 for r in touched if r >= 2 * v), excess_dims=excess)
    print('%-58s illegal gates %5d of %d %s | x %d, y %d, scratch %d on them | excess dims %d' % (
        label, len(bad), len(seq), byph, res['x_on_illegal'], res['y_on_illegal'], res['scratch_on_illegal'], excess))
    return res

out = []
# calibration: one-sided
st1 = lambda r: red([U[r]]) if r < v else ()
ce1 = lambda r: WIN if r < v else (meet(WIN, perp_fast(red([U[r - v]]))) if r < 2 * v else WIN)
out.append(run(st1, ce1, False, 'calibration: one-sided, copy per total'))
def two(Vl, label, per_reader):
    st = lambda r: red([U[r]]) if r < v else (red([Vl[r - v]]) if r < 2 * v else ())
    ce = lambda r: perp_fast(red([Vl[r]])) if r < v else (perp_fast(red([U[r - v]])) if r < 2 * v else red([1 << i for i in range(10)]))
    return run(st, ce, per_reader, label)
Vm = [S | (J9 ^ u) for u in U]
out.append(two(Vm, 'natural mirror, one copy per total', False))
out.append(two(Vm, 'natural mirror, one copy per (total, reader)', True))
tw = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'twist.json')))
if tw.get('twisted', {}).get('best'):
    A = tw['twisted']['best']['A']
    ap = lambda x: [y for y in [0]] and __import__('functools').reduce(lambda acc, i: acc ^ (A[i] if (x >> i) & 1 else 0), range(9), 0)
    Vt = [S | (J9 ^ ap(u)) for u in U]
    assert all(dot(U[b], Vt[b]) == 0 for b in range(v))
    out.append(two(Vt, 'best symmetric twisted mirror, one copy per (total, reader)', True))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'screen.json'), 'w'), indent=1)
