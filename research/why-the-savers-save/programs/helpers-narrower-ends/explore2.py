import sys; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from common import *
exec(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'explore1.py')).read().split('# min-join')[0].replace('print(', 'pass; (lambda *a, **k: None)('))
CL = c['cleanup']
# dependency DAG
writes_at = defaultdict(list); reads_at = defaultdict(list); scales_at = defaultdict(list)
pred = [set() for _ in CL]
for j, e in enumerate(CL):
    if e[0] == 'add':
        a, b = e[1], e[2]
        P = set(writes_at[b]) | set(scales_at[b]) | set(reads_at[a]) | set(scales_at[a])
        pred[j] = P
        writes_at[a].append(j); reads_at[b].append(j)
    else:
        a = e[1]
        P = set(writes_at[a]) | set(reads_at[a])
        pred[j] = P
        scales_at[a].append(j)
print('edges', sum(len(p) for p in pred))
def closure(S):
    seen = set(S); st = list(S)
    while st:
        j = st.pop()
        for i in pred[j]:
            if i not in seen: seen.add(i); st.append(i)
    return sorted(seen)
res = {}
for h in H:
    if dim(fin[h]) < 9 or not cg[h]: continue
    S = closure(cg[h])
    cur2 = {}
    def f(r): return cur2.get(r, pre[r])
    for j in S:
        e = CL[j]
        if e[0] == 'add':
            jn = join(f(e[1]), f(e[2])); cur2[e[1]] = jn; cur2[e[2]] = jn
    res[h] = (len(cg[h]), len(S), dim(pre[h]), dim(f(h)))
print('helpers considered', len(res))
print('end dim of h when its cleanup closure runs first at joins:', Counter(x[3] for x in res.values()))
print('closure sizes (min, median, max):', min(x[1] for x in res.values()), sorted(x[1] for x in res.values())[len(res)//2], max(x[1] for x in res.values()))
print('own gate counts min', min(x[0] for x in res.values()))
# partner frames directly (no closure): join of pre[h] with pre-cleanup frames of all its partners
direct = Counter()
for h in res:
    U = pre[h]
    for j in cg[h]:
        e = CL[j]
        if e[0] == 'add': U = join(U, join(pre[e[1]], pre[e[2]]))
    direct[dim(U)] += 1
print('join of h with partners pre-cleanup frames (lower bound, any order):', direct)
