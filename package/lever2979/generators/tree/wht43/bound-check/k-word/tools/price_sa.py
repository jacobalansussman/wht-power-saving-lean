import sys, os, json, gzip
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__))
from route_price import price, widths_from, IDLE
from tight_pack import cyclic_rows, validate
from sa_pack import sa
def build(profile, widths, K, v=112, copies=True, extra=None, count=72, m=73, seed=0, iters=400000):
    prof=Counter({int(k):n for k,n in profile.items()})
    w=list(widths)
    if copies: w+=[9]*v; prof[1]+=v; prof[8]+=v
    if extra:
        for r,n in extra.items(): prof[int(r)]+=n
    r=sa(w,K,m,iters=iters,seed=seed)
    if r is None: return None,None
    cost,groups,pads1=r
    rows,pads=cyclic_rows(w,groups,count,m)
    validate(rows,pads,w,count,m)
    hist=Counter({r:count*n for r,n in prof.items() if n})
    for rr,n in IDLE.items(): hist[rr]+=n*v
    hist.update(p for p in pads if p)
    W=64*v+len(pads); D=m*W-sum(r*n for r,n in hist.items())
    return dict(hist=dict(sorted(hist.items())),W=W,D=D,banks=len(pads),pads_per_replicate=sorted(Counter(pads1).items()),
                pad_cost=cost), dict(rows=rows,pads=pads,widths=w)
if __name__=='__main__':
    spec=json.loads(open(sys.argv[1]).read()); outdir=sys.argv[2]; Ks=[int(k) for k in sys.argv[3].split(',')]
    iters=int(sys.argv[4]) if len(sys.argv)>4 else 400000
    res={}
    for name,c in spec.items():
        best=None
        for K in Ks:
            for seed in range(2):
                u,sched=build(c['profile'],widths_from(c['widths']),K,copies=c.get('copies',True),extra=c.get('extra'),seed=seed,iters=iters)
                if u is None: continue
                assert u['D']==2560
                p=price(u)
                print(f"  {name} K={K} seed={seed} banks={u['banks']} W={u['W']} padcost={u['pad_cost']:.2f} pads={u['pads_per_replicate']} root={p['root_lower']} kernel={p['kernel_A']}",flush=True)
                if best is None or p['kernel_A']>best[1]['kernel_A']: best=(u,p,sched)
        u,p,sched=best
        res[name]=dict(u,**p)
        fn=os.path.join(outdir,name.replace('/','_').replace('(','_').replace(')','')+'.schedule.json.gz')
        open(fn,'wb').write(gzip.compress(json.dumps(sched,separators=(',',':')).encode(),mtime=0))
        print(f"BEST {name:40s} W={u['W']} banks={u['banks']} kernel={p['kernel_A']} dft={p['dft']}",flush=True)
    open(os.path.join(outdir,'sa.json'),'w').write(json.dumps(res,indent=1))
