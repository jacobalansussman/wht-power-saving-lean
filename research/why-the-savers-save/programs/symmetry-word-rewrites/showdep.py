import sys, json; sys.path.insert(0,'.')
import fw, search, reorder2
from ilp import bounds_PQ
W = search.Word(fw.load('K8-FULL-WORD.json.gz')); O = reorder2.OW2(W)
d = json.load(open('lbdep.json')); Pd = {int(k): tuple(v) for k, v in d['Pd'].items()}; Qd = {int(k): tuple(v) for k, v in d['Qd'].items()}
P, Q = bounds_PQ(W)
v = W.v
def role(r): return 'S%d'%r if r < v else ('T%d'%(r-v) if r < 2*v else 'H%d'%r)
r = int(sys.argv[1])
print(role(r), 'start', W.start[r], 'end dim', len(W.end[r]))
for k in O.chain[r]:
    s,i,e = W.seq[k]
    oth = [q for q in O.regs[k] if q != r]
    print(f'{s}{i:5d} {e[0]:5s} {role(e[1]):5s} {role(e[2]) if e[0]=="add" else "":6s} f{len(O.f[k])} P{len(P[k])} Q{len(Q[k])} | Pd{len(Pd[k])} Qd{len(Qd[k])}  f={O.f[k]}')
