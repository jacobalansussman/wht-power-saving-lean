"""d-halfread: screen every (second copy h3, reader c) pair of the U8 FULL-WORD for a literal half-read re-route:
c takes alpha from the merge holder a at its half line L_a before the merge (needs: c's frame just before can sit
inside L_a, c's next frame contains L_a, c idle at its old frame between), and beta*k from the partner b after the
merge at c's old read frame G (needs: b's frames after that read all contain G). Usage: python3 -I screen.py WORD copies.json"""
import sys, json
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, __import__('os').path.dirname(__file__))
from wordlib import *
w = load(sys.argv[1]); v, R = w['report']['v'], w['report']['R']; n = 2*v+R
prog = program(w); chain, bad = frame_replay(w, prog); assert not bad
C = json.load(open(sys.argv[2]))
start, end = ends(w)
cls = lambda q: 'x' if q < v else ('y' if q < 2*v else ('h-enter' if start[q] else 'h-cold'))
gates_of = defaultdict(list)
for T, (ph, e) in enumerate(prog):
    if e[0] == 'add' and ph != 'compensation':
        for q in regs(e): gates_of[q].append(T)
def frame_at(q, T):
    U = chain[q][0][1]
    for t, X in chain[q]:
        if t <= T: U = X
    return U
def next_climb(q, T):
    for t, X in chain[q]:
        if t > T: return t, X
    return None, None
def sub(A, B): return contained(A, B)
res = Counter(); per_copy = []; why_copy = Counter()
for r in C:
    Tm, phm, a, b, k = r['merge']; k = F(k); Tc, phc, h3, q = r['copy']
    P = tuple(r['P']); La = tuple(r['La']); Lb = tuple(r['Lb'])
    # a's alpha-ready time: last write into a before Tm
    Ta = max([T for T in gates_of[a] if T < Tm and prog[T][1][1] == a] or [-1])
    readers = defaultdict(list)
    for T, ph, c, cl, qq, d in r['reads']:
        if T > Tc and T < r['Tend']: readers[c].append((T, F(qq)))
    verdicts = {}
    for c, rl in readers.items():
        Tr = rl[0][0]; G = basis(prog[Tr][1][-1]) if prog[Tr][1][-1] is not None else None
        F0 = frame_at(c, Tm-1); T1, F1 = next_climb(c, Tm-1)
        lastF0 = max([T for T in gates_of[c] if T < (T1 or 10**9)] or [-1])
        if not sub(F0, La): r1 = 'R1: frame %d before the merge not inside L_a' % len(F0)
        elif F1 is None or not sub(La, F1): r1 = 'R1: next frame does not contain L_a'
        elif lastF0 >= Tm: r1 = 'R1: still busy at its frame below after the merge'
        elif max(Ta, lastF0) >= Tm: r1 = 'R1: no time window'
        else: r1 = 'ok'
        # R3: b's frames after Tr
        bad3 = [X for t, X in chain[b] if t > Tr and G is not None and not sub(G, X)]
        r3 = 'ok' if G is not None and not bad3 else 'R3: partner b climbs to %s after the read, not containing G (rank %d)' % (
            [len(X) for X in bad3][:3], len(G) if G else -1)
        verdicts[c] = (cls(c), r1, r3)
        res[(cls(c), r1, r3)] += 1
    allok = verdicts and all(x[1] == 'ok' and x[2] == 'ok' for x in verdicts.values())
    why_copy['all readers pass' if allok else 'some reader fails'] += 1
    per_copy.append(dict(P=list(P), h3=h3, a=a, b=b, verdicts={str(c): x for c, x in verdicts.items()}))
for x, y in sorted(res.items(), key=lambda z: -z[1]): print(y, x)
print(why_copy)
json.dump(per_copy, open('screen.json', 'w'), indent=0)
