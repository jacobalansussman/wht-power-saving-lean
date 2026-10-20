"""Map an annealed width-multiset schedule onto the actual helper order of Shea's C588 certificate
(schedRes in Work/GCert/Data/Gen/C588Sched.lean), rebuild cyclic rows, validate, and write a price JSON +
schedule in bank_reprice.py's format for generate_global_lean.py."""
import sys, os, re, json, gzip
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from tight_pack import cyclic_rows, validate
from route_price import IDLE
from exact_moment import certify
sched_lean, sa_sched, out = sys.argv[1], sys.argv[2], sys.argv[3]
src=open(sched_lean).read()
res=[int(x) for x in re.search(r'schedRes : List Nat := \[([0-9, ]*)\]',src).group(1).split(',')]
assert len(res)==616
s=json.loads(gzip.decompress(open(sa_sched,'rb').read()))
w=s['widths']; assert Counter(w)==Counter(res)
# recover groups (bin -> helper ids) from row 0 of the annealed schedule (bank = b*72 + (t+j)%72)
nb=len(s['pads'])//72
groups=defaultdict(list)
for q,b in enumerate(s['rows'][0]): groups[b//72].append(q)
pools=defaultdict(list)
for q,x in enumerate(res): pools[x].append(q)
G=[]
for b in range(nb):
    G.append([pools[w[q]].pop() for q in sorted(groups[b])])
assert all(not v for v in pools.values())
rows,pads=cyclic_rows(res,G,72,73)
validate(rows,pads,res,72,73)
local={1:2596,2:492,3:280,4:232,5:56,6:56,8:120}   # Shea's charged local histogram (gen-C588.json)
hist=Counter({r:72*n for r,n in local.items()})
for r,n in IDLE.items(): hist[r]+=n*112
hist.update(p for p in pads if p)
W=64*112+len(pads); mass=sum(r*n for r,n in hist.items()); assert 73*W-mass==2560
sp=out+'.schedule.json.gz'
open(sp,'wb').write(gzip.compress(json.dumps(dict(rows=rows,pads=pads,widths=res),separators=(',',':')).encode(),mtime=0))
rate=certify(dict(hist),73,W)
json.dump(dict(status='EXACT_ARITHMETIC_CANDIDATE_REPACKED_SCHEDULE',h=9,m=73,v=112,local_profile={str(k):v for k,v in local.items()},
               local_residual_count=616,banks=len(pads),W=W,deficit=2560,profile={str(k):v for k,v in sorted(hist.items())},
               pad_hist={str(k):v for k,v in sorted(Counter(pads).items())},exact_rate=rate),open(out,'w'),indent=2)
print('W',W,'banks',len(pads),'pads',sorted(Counter(pads).items()),'root',rate['root_lower_decimal'][:16])
