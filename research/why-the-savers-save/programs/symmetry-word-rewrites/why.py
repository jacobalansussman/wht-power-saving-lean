import sys; sys.path.insert(0,'.')
import fw, search
from search import le, join, meet
from collections import Counter
W = search.Word(fw.load('K8-FULL-WORD.json.gz'))
reasons = Counter()
def ok_at(k, U):
    for r in W.regs[k]:
        role = 'src' if r < W.v else ('tgt' if r < 2*W.v else 'hlp')
        if not le(U, W.end[r]): return 'end-'+role
        if not le(W.start[r], U): return 'start-'+role
        if W.sec[k] != 'T' and not le(U, W.ub[r]): return 'cone-'+role
        if k < W.nA and r in W.cframe and not le(U, W.cframe[r]): return 'centre'
    if k in W.frozen and U != W.f[k]: return 'frozen'
    return None
W.ok_at = lambda k,U: ok_at(k,U) is None
for r in range(W.n):
    role = 'src' if r < W.v else ('tgt' if r < 2*W.v else 'hlp')
    rs = W.runs(r)
    seqf = [W.start[r]] + [U for U,_ in rs] + [W.end[r]]
    for i,(U,ks) in enumerate(rs):
        for dirn, X in (('raise', seqf[i+2]), ('lower', seqf[i])):
            if X == U: continue
            # replay propagation to find first failing reason
            ch = {}; cur = lambda k: ch.get(k, W.f[k]); st = [(('up' if dirn=='raise' else 'dn'), k, X) for k in ks]; why = None
            while st and why is None:
                kind,k,Y = st.pop(); old = cur(k); new = join(old,Y) if kind=='up' else meet(old,Y)
                if new == old: continue
                why = ok_at(k,new)
                if why: break
                ch[k] = new
                for q in W.regs[k]:
                    ps = W.pos[q]; j = W.idx[(q,k)]
                    if j+1 < len(ps) and not le(new, cur(ps[j+1])): st.append(('up', ps[j+1], new))
                    if j > 0 and not le(cur(ps[j-1]), new): st.append(('dn', ps[j-1], new))
            reasons[(role, dirn, why, len(ch) > 0)] += 1
for k,c in sorted(reasons.items(), key=lambda x:-x[1]): print(k,c)
