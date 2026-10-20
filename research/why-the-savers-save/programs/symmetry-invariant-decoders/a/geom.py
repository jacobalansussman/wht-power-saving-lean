"""geom.py (t-km, part a): the face planes of a Hamming block are forced by label geometry, and the block's symmetry.
For every block b (u8_blocks.json, 14 x 8 sources) with span V_b (rank 4) and every target t that reads block b
(a0_demand.json): known(t, b) = {s in b : u_s . u_t = 0 (labels)} and u_t^perp cap V_b.  Checks (1) the demand's
support on b equals known(t, b); (2) u_t^perp cap V_b = span(known) (rank 3, so known is an affine plane); (3) the 12
foreign faces are the 12 affine planes of AG(3,2) outside one parallel class, own targets read the other class;
(4) the group of affine maps of the block (AGL(3,2) on the 8 labels, as linear maps of V_b) preserving the 12 face
needs with their sign patterns up to scale and sign changes of sources (signed permutations): its order and orbits.
usage: python -I geom.py"""
import json, itertools
from fractions import Fraction as Fr
from collections import Counter, defaultdict
d = json.load(open('a0_demand.json')); ports = d['ports']; v = len(ports)
dem = {int(t): {int(s): Fr(x) for s, x in dd.items()} for t, dd in d['demand'].items()}
blocks = json.load(open('u8_blocks.json'))['blocks']
dot = lambda a, b: bin(a & b).count('1') & 1
def basis(xs):
    b = []
    for x in xs:
        for y in b: x = min(x, x ^ y)
        if x: b.append(x); b.sort(reverse=True)
    return b
def spanset(xs):
    s = {0}
    for x in basis(xs): s |= {y ^ x for y in s}
    return s
res = Counter(); faces_per_block = []
for bi, B in enumerate(blocks):
    Vb = spanset([ports[s] for s in B]); assert len(Vb) == 16
    faces = defaultdict(set); own = defaultdict(set)
    for t in range(v):
        m = {s: x for s, x in dem[t].items() if s in B}
        if not m: res['targets not reading block'] += 1; continue
        known = {s for s in B if dot(ports[s], ports[t]) == 0}
        res['support == known'] += set(m) == known
        perp = {x for x in Vb if dot(x, ports[t]) == 0}
        res['perp cap V == span(known)'] += perp == spanset([ports[s] for s in known])
        res['rank of perp cap V = %d' % len(basis(list(perp)))] += 1
        ss = tuple(sorted(m)); x0 = m[ss[0]]
        key = (ss, tuple(m[s] / x0 for s in ss))
        (own if t in B else faces)[ss].add(key[1])
        res['own' if t in B else 'foreign'] += 1
    faces_per_block.append((len(faces), sorted(Counter(len(x) for x in faces.values()).items()), len(own)))
    # affine-plane classes: the 4-sets in the block whose labels sum to 0
    planes = [q for q in itertools.combinations(B, 4) if ports[q[0]] ^ ports[q[1]] ^ ports[q[2]] ^ ports[q[3]] == 0]
    res['affine planes per block = %d' % len(planes)] += 1
    res['faces+own planes = all affine planes'] += set(faces) | set(own) == set(planes)
print(dict(res))
print('per block (#foreign faces, patterns per face, #own planes):', Counter(map(str, faces_per_block)))
# symmetry of block 0's face-need configuration: signed affine permutations of its 8 sources
B = blocks[0]; L = {s: ports[s] for s in B}
needs = defaultdict(set)
for t in range(v):
    if t in B: continue
    m = {s: x for s, x in dem[t].items() if s in B}
    if m: needs[frozenset(m)].add(tuple(sorted((s, m[s] / m[min(m)]) for s in m)))
idx = {s: i for i, s in enumerate(B)}
pts = [L[s] for s in B]
grp = []
for pm in itertools.permutations(range(8)):
    # affine: the map pts[i] -> pts[pm[i]] must preserve the 4-sets summing to 0 (affine planes)
    ok = True
    for q in itertools.combinations(range(8), 4):
        if (pts[q[0]] ^ pts[q[1]] ^ pts[q[2]] ^ pts[q[3]] == 0) != (pts[pm[q[0]]] ^ pts[pm[q[1]]] ^ pts[pm[q[2]]] ^ pts[pm[q[3]]] == 0):
            ok = False; break
    if not ok: continue
    grp.append(pm)
print('affine permutations of the block (AGL(3,2) expected 1344):', len(grp))
faceset = {frozenset(idx[s] for s in f) for f in needs}
stabF = [pm for pm in grp if {frozenset(pm[i] for i in f) for f in faceset} == faceset]
print('stabiliser of the 12-face set:', len(stabF))
# with signs: need a sign vector eps on sources so that every face pattern maps to a face pattern (up to scale)
pat = {}
for f, ps in needs.items():
    assert len(ps) == 1
    (p_,) = ps; pat[frozenset(idx[s] for s in f)] = {idx[s]: c for s, c in p_}
def maps(pm, eps):
    for f, pp in pat.items():
        g = frozenset(pm[i] for i in f); q = pat[g]
        img = {pm[i]: pp[i] * eps[i] for i in f}
        i0 = min(img); r = img[i0] / q[i0]
        if any(img[j] != r * q[j] for j in img): return False
    return True
signed = [(pm, eps) for pm in stabF for eps in itertools.product((1, -1), repeat=8) if eps[0] == 1 and maps(pm, eps)]
print('signed symmetries of the 12 face needs (mod global sign):', len(signed))
G0perm = sorted({pm for pm, e in signed})
print('their permutation part:', len(G0perm), ' orbits on faces:',
      Counter(len({frozenset(pm[i] for i in f) for pm in G0perm}) for f in faceset))
json.dump(dict(checks=dict(res), faces_per_block=faces_per_block, agl=len(grp), stab_faces=len(stabF),
               signed_symmetries=len(signed), perm_part=len(G0perm)), open('geom.json', 'w'), indent=1)
