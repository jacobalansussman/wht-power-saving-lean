"""Price every distinct hub canonical word through Shea's route: as-is (exit-frame what-if if ends are not full),
and with non-full ends paid up to the top frame (extra climb of rank 9-d per such helper)."""
import sys, json, gzip, os
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__))
from route_price import unit, price
def read(p):
    raw=open(p,'rb').read()
    return json.loads(gzip.decompress(raw) if p.endswith('.gz') else raw)
seen={}
for p in sys.argv[2:]:
    c=read(p)
    v,h,R=c['v'],c['h'],c['R']
    if (v,h)!=(112,9): continue
    F=c['frames']
    st=[len(F[c['start'][r]]) for r in range(2*v,2*v+R)]
    fi=[len(F[c['final'][r]]) for r in range(2*v,2*v+R)]
    cst=sum(len(F[f]) for f in c['center_frames'])
    if cst!=64: continue
    prof={int(k):n for k,n in c['local_profile'].items()}
    key=json.dumps([sorted(prof.items()),sorted(Counter(zip(st,fi)).items())])
    if key in seen: seen[key]['paths'].append(p); continue
    seen[key]=dict(prof=prof,st=st,fi=fi,paths=[p],R=R)
rows=[]
for key,d in seen.items():
    prof,st,fi=d['prof'],d['st'],d['fi']
    w=[f-s for s,f in zip(st,fi)]
    full=all(f==9 for f in fi)
    u=unit(prof,w); pr=price(u,kernel=False)
    row=dict(path=d['paths'][0].split('/research/')[1],n_paths=len(d['paths']),R=d['R'],profile=prof,
             final_dims=dict(Counter(fi)),full_ended=full,width=sum(w),asis_root=pr['root_lower'],W=u['W'])
    if not full:
        extra=Counter(9-f for f in fi if f<9)
        w2=[9-s for s in st]
        u2=unit(prof,w2,extra=extra); pr2=price(u2,kernel=False)
        row.update(paid_root=pr2['root_lower'],paid_W=u2['W'])
    rows.append(row)
rows.sort(key=lambda r:-float(r['asis_root']))
json.dump(rows,open(sys.argv[1],'w'),indent=1)
print('distinct words',len(rows))
for r in rows[:25]:
    print(f"{r['asis_root']} {'FULL' if r['full_ended'] else 'exit'} paid={r.get('paid_root','-')} R={r['R']} fin={r['final_dims']} {r['path'][:90]}")
print('--- best full-ended / best paid')
best=max([r for r in rows if r['full_ended']],key=lambda r:float(r['asis_root']))
print(best)
bp=max([r for r in rows if not r['full_ended']],key=lambda r:float(r['paid_root']))
print(bp)
