import sys, json
from collections import Counter
sys.path.insert(0, '.')
import fw
w = fw.load('K8-FULL-WORD.json.gz')
p = fw.paid_profile(w); print('profile', dict(sorted(p.items())), 'cost', round(fw.cost_of_profile(p),3))
v,h,R,n,full,start,end = fw.geometry(w)
ch, start, end = fw.chains(w)
# per register chain of distinct frames
kinds = Counter(); nclimb = Counter(); ex = {}
for q in range(n):
    fr = [start[q]] + [U for k,U in ch.get(q,[])] + [end[q]]
    dims = []
    last = None
    for U in fr:
        if U != last: dims.append(len(U)); last = U
    cl = tuple(b-a for a,b in zip(dims, dims[1:]))
    role = 'src' if q < v else ('tgt' if q < 2*v else 'hlp')
    kinds[(role, cl)] += 1
    nclimb[(role, len(cl))] += 1
for (role, cl), c in sorted(kinds.items(), key=lambda x: -x[1])[:60]:
    print(role, cl, c)
print(sorted(nclimb.items()))
