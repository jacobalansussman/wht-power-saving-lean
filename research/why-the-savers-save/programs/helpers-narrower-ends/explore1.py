import sys; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from common import *
c = load(); v = c['v']; n = 2*v + c['R']
Fr = [rref(f) for f in c['frames']]
cur = [Fr[s] for s in c['start']]
runs = defaultdict(list)
def climb(r, f, ph):
    if cur[r] == f: return
    if dim(f) <= dim(cur[r]) or not inside(cur[r], f): return
    runs[r].append((dim(cur[r]), dim(f), ph)); cur[r] = f
for ph in PHASES[:-1]:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, Fr[f], ph)
        continue
    for e in c[ph]:
        if e[-1] is None: continue
        climb(e[1], Fr[e[-1]], ph)
        if e[0] == 'add': climb(e[2], Fr[e[-1]], ph)
pre = list(cur)
H = range(2*v, n)
print('dims at cleanup start: x', Counter(dim(pre[r]) for r in range(v)), 'y', Counter(dim(pre[r]) for r in range(v, 2*v)), 'h', Counter(dim(pre[r]) for r in H))
fin = [Fr[f] for f in c['final']]
print('helpers already at final at cleanup start', sum(1 for r in H if pre[r] == fin[r]), 'by final dim', Counter(dim(fin[r]) for r in H if pre[r] == fin[r]))
# cleanup gates per register
cg = defaultdict(list)
for i, e in enumerate(c['cleanup']):
    cg[e[1]].append(i)
    if e[0] == 'add': cg[e[2]].append(i)
print('helpers with no cleanup gate', sum(1 for r in H if not cg[r]))
print('cleanup gate count per helper (helpers with >=1)', Counter(len(cg[r]) for r in H if cg[r]).most_common(15))
# min-join simulation in the original order
cur2 = list(pre)
for e in c['cleanup']:
    if e[0] == 'add':
        j = join(cur2[e[1]], cur2[e[2]]); cur2[e[1]] = j; cur2[e[2]] = j
print('min-join cleanup: helper end dims', Counter(dim(cur2[r]) for r in H), ' x', Counter(dim(cur2[r]) for r in range(v)))
print('helpers ending below FULL under min-join but final FULL in U8:', sum(1 for r in H if dim(cur2[r]) < 9 and dim(fin[r]) == 9))
# 8->9 cleanup climbs (the 272)
last = {r: runs[r][-1] if runs[r] else None for r in H}
