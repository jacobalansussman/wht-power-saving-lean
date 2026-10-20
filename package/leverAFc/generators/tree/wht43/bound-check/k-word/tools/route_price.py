"""Price a local word through Shea's complex588 route (source copies + #401/#377 dedicated banks,
cyclic schedule, finite fill s=40, Sussman's rate generator), from (local profile, helper widths).
Functions copied from shea256/fourier-transform-below-nlogn (bank_reprice.py, generate_global_lean.py);
foldlib/gen_rates from the same tree (Sussman's rate generator). Arithmetic only."""
import sys, json
from collections import Counter
from fractions import Fraction as Q
sys.path.insert(0, __import__('os').path.dirname(__file__))
from exact_moment import certify
import foldlib as fl

IDLE = {48: 32, 32: 32, 16: 32, 8: 48}   # GX.bankPrice377 per source class (Bridge377/Bank.lean)

def cyclic_schedule(widths, count=72, m=73):
    remaining=sorted(range(len(widths)),key=lambda q:(-widths[q],q))
    groups=[]
    while remaining:
        dp={0:[]}
        for q in remaining:
            for load in sorted(list(dp),reverse=True):
                s=load+widths[q]
                if s<=m and s not in dp: dp[s]=dp[load]+[q]
            if m in dp: break
        best=dp[max(dp)]
        assert best and len(best)<=count
        groups.append(best)
        taken=set(best)
        remaining=[q for q in remaining if q not in taken]
    pads=[]
    for b,group in enumerate(groups):
        pads.extend([m-sum(widths[q] for q in group)]*count)
    return groups,pads

def unit(profile, widths, v=112, cst=64, copies=True, idle=IDLE, calls=72, m=73, extra=None):
    prof=Counter({int(k):n for k,n in profile.items()})
    w=list(widths)
    if copies:
        w += [9]*v
        prof[1]+=v; prof[8]+=v
    if extra:
        for r,n in extra.items(): prof[int(r)]+=n
    groups,pads=cyclic_schedule(w,calls,m)
    hist=Counter({r:calls*n for r,n in prof.items() if n})
    for r,n in idle.items(): hist[r]+=n*v
    hist.update(p for p in pads if p)
    W=64*v+len(pads)
    mass=sum(r*n for r,n in hist.items())
    D=m*W-mass
    return dict(hist=dict(sorted(hist.items())),W=W,D=D,banks=len(pads),bins=len(groups),
                pad_hist=dict(sorted(Counter(pads).items())),local_mass=sum(r*n for r,n in prof.items()),
                width=sum(w))

def price(u, m=73, kernel=True):
    rate=certify(u['hist'],m,u['W'])
    out=dict(root_lower=rate['root_lower_decimal'][:16])
    if kernel:
        uu=dict(m=m,W=u['W'],D=u['D'],H=u['hist'])
        A=int(Q(rate['root_lower'])*10**10)
        while fl.plan_block(uu,A,10,40) is None: A-=1
        out.update(kernel_A=A, dft=A-1)
    return out

def widths_from(d):
    w=[]
    for k,n in sorted(((int(k),n) for k,n in d.items()), reverse=True): w += [k]*n
    return w

if __name__=='__main__':
    spec=json.loads(open(sys.argv[1]).read())
    res={}
    for name,c in spec.items():
        u=unit(c['profile'],widths_from(c['widths']),copies=c.get('copies',True),extra=c.get('extra'))
        if u['D']!=2560: print(name,'WARNING deficit',u['D'])
        p=price(u,kernel=c.get('kernel',True))
        res[name]=dict(u,**p)
        print(f"{name:40s} W={u['W']} banks={u['banks']} width={u['width']} D={u['D']} root={p['root_lower']} kernel={p.get('kernel_A')} dft={p.get('dft')}")
    if len(sys.argv)>2: open(sys.argv[2],'w').write(json.dumps(res,indent=1))
