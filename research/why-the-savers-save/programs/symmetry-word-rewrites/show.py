import sys; sys.path.insert(0,'.')
import fw, search
W = search.Word(fw.load('K8-FULL-WORD.json.gz'))
v = W.v
def role(r): return 'S%d'%r if r < v else ('T%d'%(r-v) if r < 2*v else 'H%d'%r)
def fs(U): return '<'+','.join(str(x) for x in U)+'>'
r = int(sys.argv[1])
print(role(r), 'start', fs(W.start[r]), 'end', fs(W.end[r]), 'ub', len(W.ub[r]))
last = W.start[r]
for k in W.pos[r]:
    s,i,e = W.seq[k]
    other = [q for q in W.regs[k] if q != r]
    U = W.f[k]
    mark = ('climb+%d'%(len(U)-len(last))) if U != last else ''
    print(f'{s}{i:5d} {e[0]:5s} {role(e[1]):5s} {"<- "+role(e[2]) if e[0]=="add" else "":10s} q={e[3]}/{e[4]} dim{len(U)} {mark:9s} {fs(U)}')
    last = U
