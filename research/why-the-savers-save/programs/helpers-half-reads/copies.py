"""d-halfread: census of the 168 second copies (h3 += q*a at P after the merge a += k*b at P) of the U8 FULL-WORD:
writers into h3, readers of h3 (who, phase, frame rank, inside the packet window), consumers of the packet.
Usage: python3 -I copies.py FULL-WORD.json.gz out.json"""
import sys, json
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, __import__('os').path.dirname(__file__))
from wordlib import *
w = load(sys.argv[1]); v, R = w['report']['v'], w['report']['R']; n = 2*v+R
prog = program(w); chain, bad = frame_replay(w, prog); assert not bad
cls = lambda q: 'x' if q < v else ('y' if q < 2*v else 'h')
rows = [{k: F(1)} for k in range(n)]; readval = {}
for T, (ph, e) in enumerate(prog):
    if e[0] == 'add': readval[T] = dict(rows[e[2]])
    apply(e, rows)
def frame_at(q, T):
    U = chain[q][0][1]
    for t, X in chain[q]:
        if t <= T: U = X
    return U
planes = defaultdict(list)
for q in range(n):
    for T, U in chain[q]:
        if len(U) == 3: planes[U].append((q, T))
start, end = ends(w)
out = []; stat = Counter()
for P in sorted(planes):
    gP = [(T, ph, e) for T, (ph, e) in enumerate(prog) if e[0] == 'add' and ph in ('A', 'B', 'tail') and basis(e[-1]) == P]
    (Tm, phm, em), (Tc, phc, ec) = gP
    a, b, k = em[1], em[2], F(em[3], em[4]); h3, q = ec[1], F(ec[3], ec[4]); assert ec[2] == a
    La, Lb = frame_at(a, Tm-1), frame_at(b, Tm-1)
    alpha = {s: x for s, x in readval[Tc].items()}  # a's value at the copy (= p + a's dirt)
    ev = [(T, ph, e) for T, (ph, e) in enumerate(prog) if h3 in regs(e) and e[0] == 'add' and ph != 'compensation']
    comp = sum(1 for ph, e in prog if ph == 'compensation' and h3 in regs(e))
    writes = [(T, ph, e[2], str(F(e[3], e[4])), len(basis(e[-1])) if e[-1] is not None else None) for T, ph, e in ev if e[1] == h3]
    reads = [(T, ph, e[1], cls(e[1]), str(F(e[3], e[4])), len(basis(e[-1])) if e[-1] is not None else None) for T, ph, e in ev if e[2] == h3]
    # packet window: from Tc to the write that removes a's packet (first later write from a with -q... or any write from a)
    Tend = min([T for T, ph, s, *_ in writes if T > Tc and s == a] or [10**9])
    inwin = [r for r in reads if Tc < r[0] < Tend]
    other_w = [x for x in writes if x[0] != Tc and x[0] != Tend]
    # net packet coefficient a consumer gets in the window
    cons = Counter()
    for T, ph, c, cl, qq, d in inwin: cons[c] += F(qq)
    pre = Counter()
    for T, ph, c, cl, qq, d in reads:
        if T < Tc: pre[c] += F(qq)
    rec = dict(P=list(P), merge=[Tm, phm, a, b, str(k)], copy=[Tc, phc, h3, str(q)], La=list(La), Lb=list(Lb),
               h3_start=list(start[h3]), h3_end=len(end[h3]), comp_reads=comp, writes=writes, reads=reads, Tend=Tend,
               window_reads=inwin, other_writes=other_w,
               consumers={str(c): str(x) for c, x in cons.items() if x}, pre_reads={str(c): str(x) for c, x in pre.items() if x})
    out.append(rec)
    key = ('entering' if start[h3] else 'cold', 'other_writes=%d' % len(other_w),
           'consumers ' + ','.join(sorted(cls(int(c)) for c in rec['consumers'])))
    stat[key] += 1
for k_, x in sorted(stat.items()): print(x, k_)
json.dump(out, open(sys.argv[2], 'w'), indent=0)
