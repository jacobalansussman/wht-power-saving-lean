import json,sys,gzip
from collections import Counter
def read(p):
    raw=open(p,'rb').read()
    return json.loads(gzip.decompress(raw) if p.endswith('.gz') else raw)
for p in sys.argv[1:]:
    c=read(p)
    v,h,R=c['v'],c['h'],c['R']
    F=c['frames']
    fin=Counter(len(F[c['final'][r]]) for r in range(2*v,2*v+R))
    st=Counter(len(F[c['start'][r]]) for r in range(2*v,2*v+R))
    w=Counter(len(F[c['final'][r]])-len(F[c['start'][r]]) for r in range(2*v,2*v+R))
    prof={int(k):n for k,n in c['local_profile'].items()}
    print(p)
    print('  v',v,'h',h,'R',R,'profile',dict(sorted(prof.items())),'mass',sum(k*n for k,n in prof.items()))
    print('  helper final dims',dict(sorted(fin.items())),' start dims',dict(sorted(st.items())))
    print('  widths',dict(sorted(w.items())),'sum',sum(k*n for k,n in w.items()), 'cst',sum(len(F[f]) for f in c['center_frames']))
