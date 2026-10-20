import sys; sys.path.insert(0,'.')
import fw, search
from search import le, join, meet
from collections import Counter
W = search.Word(fw.load(sys.argv[1] if len(sys.argv)>1 else 'K8-FULL-WORD.json.gz'))
N = len(W.seq)
# forward lazy P
cur = list(W.start); P = [None]*N
for k in range(N):
    if k == W.nA:
        for q,U in W.cframe.items(): pass
    rr = W.regs[k]; U = ()
    for r in rr: U = join(U, cur[r])
    P[k] = U
    for r in rr: cur[r] = U
# backward Q
cur = list(W.end); Q = [None]*N
for k in reversed(range(N)):
    rr = W.regs[k]; U = W.full
    for r in rr:
        U = meet(U, cur[r])
        if W.sec[k] != 'T': U = meet(U, W.ub[r])
        if k < W.nA and r in W.cframe: U = meet(U, W.cframe[r])
    Q[k] = U
    for r in rr: cur[r] = U
stat = Counter()
for k in range(N):
    f = W.f[k]
    stat[(W.sec[k], 'P<=f' if le(P[k],f) else 'P!<=f', 'f<=Q' if le(f,Q[k]) else 'f!<=Q', 'forced' if P[k]==Q[k] else ('f=P' if f==P[k] else ('f=Q' if f==Q[k] else 'mid')))] += 1
for k,c in sorted(stat.items()): print(k,c)

def cost_with(fr):
    g = lambda k: fr[k]
    return sum(W.chain_cost(r, g) for r in range(W.n))
from collections import Counter as C
def prof(fr):
    g = lambda k: fr[k]; p = C()
    for r in range(W.n): p.update(W.chain_profile(r, g))
    return dict(sorted(p.items()))
print('orig', round(cost_with(W.f),3), prof(W.f))
print('allP', round(cost_with(P),3), prof(P))
print('allQ', round(cost_with(Q),3), prof(Q))
# centre check for P and Q
for name, fr in (('P',P),('Q',Q)):
    bad = 0
    for q,U in W.cframe.items():
        pa = [k for k in W.pos[q] if k < W.nA]
        last = fr[pa[-1]] if pa else W.start[q]
        if last != U: bad += 1
    print(name, 'centre mismatches', bad)
