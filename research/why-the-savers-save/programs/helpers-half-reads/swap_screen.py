"""d-halfread: screen for the 'copy-of-copy' re-route. Many second copies h3 (of P) are read by c, the second copy of
another plane P2, which then holds p2 + p.  Re-route: c stops visiting P2; it takes alpha of P at L_a (pre-merge),
then at its old frame G reads b (k*beta) and a2 (the holder of p2 after P2's merge).  h3 retires; P and P2 both keep
2 visitors.  Tests the frame conditions literally (c inside L_a before Tm with L_a inside G; b's and a2's later frames
contain G; a2 still at a frame inside G). Usage: python3 -I swap_screen.py WORD copies.json"""
import sys, json
from collections import Counter
sys.path.insert(0, __import__('os').path.dirname(__file__))
from wordlib import *
w = load(sys.argv[1]); v, R = w['report']['v'], w['report']['R']; n = 2*v+R
prog = program(w); chain, bad = frame_replay(w, prog); C = json.load(open(sys.argv[2]))
bycopy = {r['copy'][2]: r for r in C}
def frame_at(q, T):
    U = chain[q][0][1]
    for t, X in chain[q]:
        if t <= T: U = X
    return U
def later_ok(q, T, G):
    return all(contained(G, X) for t, X in chain[q] if t > T)
stat = Counter(); passes = []
for r in C:
    Tm, _, a, b, k = r['merge']; Tc, _, h3, q = r['copy']; La = tuple(r['La']); P = tuple(r['P'])
    rd = [(T, c) for T, ph, c, cl, qq, d in r['reads'] if Tc < T < r['Tend']]
    for T, c in rd:
        if c not in bycopy: stat['reader is not a second copy'] += 1; continue
        r2 = bycopy[c]; Tm2, _, a2, b2, k2 = r2['merge']; Tc2 = r2['copy'][0]
        G = basis(prog[T][1][-1]); conds = []
        F0 = frame_at(c, Tm-1)
        conds.append(('c inside L_a before merge', contained(F0, La)))
        conds.append(('L_a inside G', contained(La, G)))
        conds.append(('c reaches L_a before merge (c not at P2 yet)', Tc2 > Tm or True))
        conds.append(('b at a frame inside G at T', contained(frame_at(b, T), G)))
        conds.append(("b's later frames contain G", later_ok(b, T, G)))
        conds.append(('a2 at a frame inside G at T', contained(frame_at(a2, T), G)))
        conds.append(("a2's later frames contain G", later_ok(a2, T, G)))
        first = next((nm for nm, ok in conds if not ok), 'all pass')
        stat[first] += 1
        if first == 'all pass': passes.append((P, h3, c, a2))
for s, x in sorted(stat.items(), key=lambda z: -z[1]): print(x, s)
print('passes', passes[:10])
# ---- details of the passes: value checks and busy intervals
from fractions import Fraction as Fr
rows = [{kk: Fr(1)} for kk in range(n)]; valat = {}
watch = set()
for P, h3, c, a2 in passes:
    r = [x for x in C if x['copy'][2] == h3][0]; watch |= {r['merge'][2], r['merge'][3], a2, c}
snap = {q: [] for q in watch}
for T, (ph, e) in enumerate(prog):
    for q in regs(e):
        if q in watch: snap[q].append((T, dict(rows[q])))
    apply(e, rows)
    for q in regs(e):
        if q in watch: snap[q].append((T + 0.5, dict(rows[q])))
def val(q, T):  # value just before time T
    best = {q: Fr(1)}
    for t, d in snap[q]:
        if t < T: best = d
    return best
cl = lambda q: 'x' if q < v else ('y' if q < 2*v else ('h-enter' if start[q] else 'h-cold'))
start, end = ends(w)
gl = {}
for T, (ph, e) in enumerate(prog):
    if ph == 'compensation': continue
    for q in regs(e): gl.setdefault(q, []).append(T)
print('\nP h3 c a2 | a b classes | a2 same value at Tc2 and T | b holds only beta at T | c gates between (Tm, T) excluding h3/a2 | readers of c in (Tc2, T)')
uc = Counter()
for P, h3, c, a2 in passes:
    r = [x for x in C if x['copy'][2] == h3][0]; Tm, _, a, b, k = r['merge']
    r2 = bycopy[c]; Tc2 = r2['copy'][0]
    T = [t for t, ph, cc, *_ in r['reads'] if cc == c and t > r['copy'][0]][0]
    same_a2 = val(a2, Tc2) == val(a2, T)
    bv = val(b, T); srcb = {s: x for s, x in bv.items() if s < v}; junkb = {s: x for s, x in bv.items() if s >= v}
    cg = [t for t in gl[c] if Tm < t < T and prog[t][1][0] == 'add' and h3 not in regs(prog[t][1]) and a2 not in regs(prog[t][1])]
    rc = [t for t in gl[c] if Tc2 < t < T and prog[t][1][0] == 'add' and prog[t][1][2] == c]
    uc[c] += 1
    print(P, h3, c, a2, '|', cl(a), cl(b), '|', same_a2, '| b junk', len(junkb), '|', [(t, prog[t][0], prog[t][1][1], prog[t][1][2]) for t in cg][:4], '|', rc[:4])
print('distinct c', len(uc), 'distinct h3', len({p[1] for p in passes}))
# ---- strict timing: c must be idle between the pre-merge read and its read of h3 (apart from gates that go)
print('\nstrict timing screen')
st = Counter(); keep = []
for P, h3, c, a2 in passes:
    r = [x for x in C if x['copy'][2] == h3][0]; Tm, _, a, b, k = r['merge']
    r2 = bycopy[c]; Tc2 = r2['copy'][0]
    T = [t for t, ph, cc, *_ in r['reads'] if cc == c and t > r['copy'][0]][0]
    Ta = max([t for t in gl[a] if t < Tm and prog[t][1][0] == 'add' and prog[t][1][1] == a] or [-1])
    lastc = max([t for t in gl[c] if t < Tm] or [-1])
    mid = [t for t in gl[c] if Tm < t < T and h3 not in regs(prog[t][1]) and t != Tc2]
    why = 'T1 no window before the merge' if max(Ta, lastc) >= Tm - 1 else ('T2 c busy between the merge and its read of h3: %d gates' % len(mid) if mid else 'pass')
    st[why] += 1
    if why == 'pass': keep.append((P, h3, c, a2))
for s, x in st.items(): print(x, s)
print(keep)
print('\none-gate cases')
for P, h3, c, a2 in passes:
    r = [x for x in C if x['copy'][2] == h3][0]; Tm, phm, a, b, k = r['merge']
    r2 = bycopy[c]; Tc2 = r2['copy'][0]
    T = [t for t, ph, cc, *_ in r['reads'] if cc == c and t > r['copy'][0]][0]
    mid = [t for t in gl[c] if Tm < t < T and h3 not in regs(prog[t][1]) and t != Tc2]
    if len(mid) <= 2:
        print(P, 'h3', h3, 'c', c, 'merge', Tm, phm, 'a', a, 'b', b, 'Tc2', Tc2, 'T', T, 'mid', [(t, prog[t][0], prog[t][1][:5], len(basis(prog[t][1][-1])) if prog[t][1][-1] is not None else None) for t in mid])
