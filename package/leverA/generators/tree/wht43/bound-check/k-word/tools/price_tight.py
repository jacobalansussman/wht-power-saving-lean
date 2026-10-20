import sys, os, json, gzip
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__))
from route_price import price, widths_from, IDLE
from tight_pack import split_exact, groups_from, cyclic_rows, validate
from exact_moment import certify

def build(profile, widths, v=112, copies=True, extra=None, count=72, m=73, seed=1):
    prof=Counter({int(k):n for k,n in profile.items()})
    w=list(widths)
    if copies: w+=[9]*v; prof[1]+=v; prof[8]+=v
    if extra:
        for r,n in extra.items(): prof[int(r)]+=n
    sol=split_exact(w,m,seed=seed)
    assert sol, 'no exact split found'
    sizes,chosen,leftover=sol
    groups,lg=groups_from(w,sizes,chosen,leftover)
    rows,pads=cyclic_rows(w,groups+[lg],count,m)
    validate(rows,pads,w,count,m)
    hist=Counter({r:count*n for r,n in prof.items() if n})
    for r,n in IDLE.items(): hist[r]+=n*v
    hist.update(p for p in pads if p)
    W=64*v+len(pads); D=m*W-sum(r*n for r,n in hist.items())
    return dict(hist=dict(sorted(hist.items())),W=W,D=D,banks=len(pads),pad_hist=dict(Counter(pads)),
                leftover=leftover,local_profile=dict(prof)), dict(rows=rows,pads=pads,widths=w)

if __name__=='__main__':
    spec=json.loads(open(sys.argv[1]).read()); outdir=sys.argv[2]
    res={}
    for name,c in spec.items():
        u,sched=build(c['profile'],widths_from(c['widths']),copies=c.get('copies',True),extra=c.get('extra'))
        p=price(u)
        res[name]=dict(u,**p)
        fn=os.path.join(outdir,name.replace('/','_').replace('(','_').replace(')','')+'.schedule.json.gz')
        open(fn,'wb').write(gzip.compress(json.dumps(sched,separators=(',',':')).encode(),mtime=0))
        print(f"{name:40s} W={u['W']} banks={u['banks']} pads={u['pad_hist']} leftover={u['leftover']} D={u['D']} root={p['root_lower']} kernel={p['kernel_A']} dft={p['dft']}")
    open(os.path.join(outdir,'tight.json'),'w').write(json.dumps(res,indent=1))
