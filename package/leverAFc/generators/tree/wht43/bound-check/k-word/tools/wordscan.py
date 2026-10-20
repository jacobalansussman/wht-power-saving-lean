"""Read-only analysis of a hub canonical word (CANONICAL-WORD.json[.gz]): local profile, helper (entrance, exit)
dims, nesting E <= U over GF(2)^9, width = dim U - dim E, which source/target rows are written, gate counts."""
import sys, json, gzip
from collections import Counter
def read(p):
    raw=open(p,'rb').read()
    return json.loads(gzip.decompress(raw) if p.endswith('.gz') else raw)
def rank(vecs):
    basis=[]  # gf2 elimination on ints
    for v in vecs:
        for b in basis:
            v=min(v,v^b)
        if v: basis.append(v)
    return len(basis)
def contains(U,E):  # span(E) <= span(U)
    return rank(list(U)+list(E))==rank(U)
def analyse(p):
    c=read(p); v,h,R=c['v'],c['h'],c['R']; F=c['frames']
    st,fi=c['start'],c['final']
    out=dict(path=p,v=v,h=h,R=R,profile={int(k):n for k,n in c['local_profile'].items()})
    for k in F: assert rank(k)==len(k), 'frame not independent'
    pairs=[];bad=[]
    for r in range(2*v,2*v+R):
        E,U=F[st[r]],F[fi[r]]
        if not contains(U,E): bad.append(r)
        pairs.append((len(E),len(U)))
    out['helper_pairs']=sorted(Counter(pairs).items())
    out['nest_fail']=bad
    out['widths']=dict(sorted(Counter(u-e for e,u in pairs).items()))
    out['width']=sum(u-e for e,u in pairs)
    out['src_dims']=sorted(Counter((len(F[st[r]]),len(F[fi[r]])) for r in range(v)).items())
    out['tgt_dims']=sorted(Counter((len(F[st[r]]),len(F[fi[r]])) for r in range(v,2*v)).items())
    out['cst']=sum(len(F[f]) for f in c['center_frames'])
    wr={}
    for ph in ['precondition','A','cut','compensation','B','cleanup','tail']:
        g=c.get(ph) or []
        srcw=Counter(x[1] for x in g if x[0] in ('add','scale') and x[1]<v)
        tgtw=Counter(x[1] for x in g if x[0] in ('add','scale') and v<=x[1]<2*v)
        ops=Counter(x[0] for x in g)
        wr[ph]=dict(n=len(g),ops=dict(ops),src_rows_written=len(srcw),tgt_rows_written=len(tgtw))
    out['phases']=wr
    allsrc=set()
    for ph in ['precondition','A','cut','compensation','B','cleanup','tail']:
        for x in c.get(ph) or []:
            if x[0] in ('add','scale') and x[1]<v: allsrc.add(x[1])
    out['sources_mutated']=len(allsrc)
    mass=sum(k*n for k,n in out['profile'].items())
    out['mass']=mass
    out['mass_check']=out['width']+sum(f-s for (s,f),n in out['src_dims'] for _ in range(n))+sum(f-s for (s,f),n in out['tgt_dims'] for _ in range(n))+out['cst']
    return out
if __name__=='__main__':
    res=[analyse(p) for p in sys.argv[1:]]
    for o in res:
        print(json.dumps({k:o[k] for k in o if k!='phases'}))
        print('  phases',json.dumps(o['phases']))
