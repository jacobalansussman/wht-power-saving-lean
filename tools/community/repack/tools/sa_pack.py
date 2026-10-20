"""Simulated-annealing bin packing for the cyclic dedicated-bank schedule: per replicate, K bins of capacity 73,
minimise sum_b f(pad_b), f(p) = p ln(73/p) (the first-order moment cost of a padding block of rank p)."""
import math, random
from collections import Counter
def f(p,m=73): return 0.0 if p==0 else p*math.log(m/p)
def sa(widths, K, m=73, iters=400000, seed=0, T0=3.0, T1=0.01):
    rng=random.Random(seed)
    n=len(widths)
    order=sorted(range(n),key=lambda q:-widths[q])
    load=[0]*K; where=[0]*n
    for q in order:   # first fit decreasing into K bins, overflow into least loaded
        b=min(range(K),key=lambda b:(load[b]+widths[q]>m, -load[b] if load[b]+widths[q]<=m else load[b]))
        where[q]=b; load[b]+=widths[q]
    def cost_bin(l): return f(m-l) if l<=m else 1000+50*(l-m)
    cost=sum(cost_bin(l) for l in load)
    best=cost; bestw=where[:]
    for it in range(iters):
        T=T0*(T1/T0)**(it/iters)
        if rng.random()<0.5:
            q=rng.randrange(n); a=where[q]; b=rng.randrange(K)
            if a==b: continue
            w=widths[q]
            d=cost_bin(load[a]-w)+cost_bin(load[b]+w)-cost_bin(load[a])-cost_bin(load[b])
            if d<=0 or rng.random()<math.exp(-d/T):
                load[a]-=w; load[b]+=w; where[q]=b; cost+=d
        else:
            q=rng.randrange(n); r=rng.randrange(n); a=where[q]; b=where[r]
            if a==b or widths[q]==widths[r]: continue
            dw=widths[r]-widths[q]
            d=cost_bin(load[a]+dw)+cost_bin(load[b]-dw)-cost_bin(load[a])-cost_bin(load[b])
            if d<=0 or rng.random()<math.exp(-d/T):
                load[a]+=dw; load[b]-=dw; where[q]=b; where[r]=a; cost+=d
        if cost<best-1e-12: best=cost; bestw=where[:]
    load=[0]*K
    for q,b in enumerate(bestw): load[b]+=widths[q]
    if max(load)>m: return None
    groups=[[q for q in range(n) if bestw[q]==b] for b in range(K)]
    groups=[g for g in groups if g]
    return best, groups, sorted(m-l for l in load if l)
