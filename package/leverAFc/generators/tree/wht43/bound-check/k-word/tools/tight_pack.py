"""Tighter bank packing for Shea's cyclic dedicated-bank schedule.
Per call: split the helper widths into (nbins-1) groups of EXACTLY m=73 plus one leftover group, so that the
padding is ONE block of rank (slack) per replicate instead of many rank-1 blocks.  Optionally pack the 72
leftover copies across calls into the fewest banks (each bank at most one helper per call)."""
import random, sys, json, os
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__))

def patterns(sizes, m=73):
    out=[]
    def rec(i,rem,cur):
        if i==len(sizes):
            if rem==0: out.append(tuple(cur))
            return
        s=sizes[i]
        for k in range(rem//s,-1,-1):
            rec(i+1,rem-k*s,cur+[k])
    rec(0,m,[])
    return out

def split_exact(widths, m=73, trials=20000, seed=1):
    cnt=Counter(widths); sizes=sorted(cnt,reverse=True)
    total=sum(widths); nfull=total//m; left=total-nfull*m
    pats=patterns(sizes,m)
    rng=random.Random(seed)
    for t in range(trials):
        rem=[cnt[s] for s in sizes]; chosen=[]
        ok=True
        for b in range(nfull):
            feas=[p for p in pats if all(p[i]<=rem[i] for i in range(len(sizes)))]
            if not feas: ok=False; break
            # prefer patterns that use the scarcest items proportionally (random tie-break)
            def score(p):
                return sum(p[i]/(rem[i]+1e-9) for i in range(len(sizes)) if p[i])+rng.random()*0.5
            p=max(feas,key=score) if rng.random()<0.7 else rng.choice(feas)
            for i in range(len(sizes)): rem[i]-=p[i]
            chosen.append(p)
        if ok:
            leftover={sizes[i]:rem[i] for i in range(len(sizes)) if rem[i]}
            assert sum(s*n for s,n in leftover.items())==left
            return sizes,chosen,leftover
    return None

def groups_from(widths, sizes, chosen, leftover):
    """turn pattern counts into index groups (helper indices)"""
    pools={s:[q for q in range(len(widths)) if widths[q]==s] for s in sizes}
    groups=[]
    for p in chosen:
        g=[]
        for i,s in enumerate(sizes):
            for _ in range(p[i]): g.append(pools[s].pop())
        groups.append(g)
    g=[]
    for s,n in leftover.items():
        for _ in range(n): g.append(pools[s].pop())
    assert all(not v for v in pools.values())
    return groups, g

def cyclic_rows(widths, groups, count=72, m=73):
    rows=[[None]*len(widths) for _ in range(count)]; pads=[]
    for b,group in enumerate(groups):
        pads.extend([m-sum(widths[q] for q in group)]*count)
        for t in range(count):
            for j,q in enumerate(group): rows[t][q]=b*count+(t+j)%count
    return rows,pads

def validate(rows,pads,widths,count=72,m=73):
    rec=[0]*len(pads)
    assert len(rows)==count
    for row in rows:
        assert len(row)==len(widths) and len(set(row))==len(row)
        for q,b in enumerate(row):
            assert 0<=b<len(pads); rec[b]+=widths[q]
    assert all(0<=p<m and l+p==m for l,p in zip(rec,pads))
