"""d-halfread BIT SIDE: retire extra carrier 470 of #589's tuned bit word (research/e8-bit-p111-tuned/tune).
Carrier 470 copies source x77 at x77's line frame (A gate 341) only so that helper 444 can read the row x77 at frame
79 (A gate 2104); by then x77 itself sits at frame 186 = [33,77] and is unchanged until B.  Here 444 reads x77 directly
at frame 186 (444: 0 -> [33,77] -> 509 instead of 0 -> [77] -> 509; 509 contains 186), the C undo of 444 reads x77
too, carrier 470's own gates (A 341, C 5242) go, 470 leaves the dirt section D, and 470 stays an unused wire (final
frame 0, no climb).  Exact replay with the hub's wlib.replay (F2 endpoint on every column, nested frames).
Usage: python3 -I retire470.py"""
import sys, os, json, copy
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import wlib
w0 = wlib.load(os.path.join(HERE, 'bitE8mut-p111-tuned.json.gz'))
base = wlib.paid(w0)
print('base local paid profile', base, 'mass', sum(r*n for r, n in base.items()))
w = copy.deepcopy(w0)
fr = w['frames']; f186 = fr.index([33, 77]) if [33, 77] in fr else None
print('frame [33,77] index', f186)
A = w['A']; C = w['C']
i341 = [i for i, g in enumerate(A) if g == [79, 470, [77]]]; i2104 = [i for i, g in enumerate(A) if g == [79, 444, [470]]]
c_undo = [i for i, g in enumerate(C) if g[1] == 444 and g[2] == [470]]; c_470 = [i for i, g in enumerate(C) if g[1] == 470 and g[2] == [77]]
print('A gates', i341, i2104, 'C gates', c_undo, c_470)
A[i2104[0]] = [f186, 444, [77]]
C[c_undo[0]] = [C[c_undo[0]][0], 444, [77]]
for i in sorted(i341, reverse=True): del A[i]
for i in sorted(c_470, reverse=True): del C[i]
nd = 0
for g in w['D']:
    if 470 in g[2]: g[2].remove(470); nd += 1
w['D'] = [g for g in w['D'] if g[2]]
w['final'][470] = 0
print('removed 470 from', nd, 'D gates; 470 gates left:', sum(1 for sec in 'DABC' for g in w[sec] if g[1] == 470 or 470 in g[2]))
try:
    new = wlib.paid(w); ok = True
except AssertionError as e:
    ok = False; new = None; print('REPLAY FAILED at', e)
print('exact replay (endpoint on every column, nested frames):', ok)
if ok:
    print('new local paid profile', new, 'mass', sum(r*n for r, n in new.items()))
    d = {r: new.get(r, 0) - base.get(r, 0) for r in set(base) | set(new)}
    print('delta', {r: x for r, x in sorted(d.items()) if x})
    wlib.save(w, os.path.join(HERE, 'bitE8mut-p111-tuned-retire470.json.gz'))
    json.dump(dict(base=base, new=new, delta={r: x for r, x in d.items() if x}), open(os.path.join(HERE, 'retire470.json'), 'w'), indent=1)
# ---- renumbered word (register 470 deleted, R = 438) for the hub's independent checker xcheck_p.py
w2 = copy.deepcopy(w)
X = 470
f = lambda q: q - 1 if q > X else q
for sec in 'DABC':
    w2[sec] = [[g[0], f(g[1]), [f(s) for s in g[2]]] for g in w2[sec]]
w2['ret'] = [f(q) for q in w2['ret']]
del w2['start'][X]; del w2['final'][X]
w2['R'] = w0['R'] - 1
w2['stats'] = dict(w2['stats']); w2['stats']['extra_carriers'] = w2['stats']['extra_carriers'] - 1
assert w2['C'] == wlib.rebuild_C(w2), 'C is not the literal reverse'
new2 = wlib.paid(w2)
print('renumbered: R', w2['R'], 'replay ok, profile', new2, 'C is the literal reverse of A,B z-gates:', True)
sha = wlib.save(w2, os.path.join(HERE, 'bitE8mut-p111-tuned-retire470-R438.json.gz')); print('sha256', sha)
