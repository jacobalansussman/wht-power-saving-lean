"""v-s-fix/indep_local.py: independent recount of the local-reading 'knows' pattern on the #589 bit ports,
plus the hub's own projector identities re-verified exactly (Fraction), for every port.
usage: python3 -I indep_local.py  (reads ../s-fix/GEOMETRY-nine-bit589.json; sha256 checked)"""
import json, hashlib
from fractions import Fraction as F
p='../../data/hub-inputs/GEOMETRY-nine-bit589.json'
print('sha256', hashlib.sha256(open(p,'rb').read()).hexdigest())
g=json.load(open(p)); ports=g['ports']; V=len(ports)
q=[[F(x) for x in t['q']] for t in ports]
dot=lambda u,v: sum(a*b for a,b in zip(u,v))
norms={dot(x,x) for x in q}
off=sum(1 for s in range(V) for t in range(V) if s!=t and dot(q[s],q[t])==0)
deg={sum(1 for t in range(V) if t!=s and dot(q[s],q[t])==0) for s in range(V)}
print('ports',V,'norms q.q',sorted(norms),'ordered orthogonal off-diagonal pairs',off,'of',V*(V-1),'degree set',sorted(deg))
# the hub identities rho*diag(X,I7)*rho for X in P,H,0,K (exact)
def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
bad=0
for t in ports:
    R=[[F(x) for x in r] for r in t['rho15']]
    qq=[F(x) for x in t['q']]; n2=dot(qq,qq)
    P=[[qq[i]*qq[j]/n2 for j in range(8)] for i in range(8)]
    I8=[[F(int(i==j)) for j in range(8)] for i in range(8)]
    K=[[I8[i][j]-P[i][j] for j in range(8)] for i in range(8)]
    Z=[[F(0)]*8 for _ in range(8)]
    def diag(A,b):
        M=[[F(0)]*15 for _ in range(15)]
        for i in range(8):
            for j in range(8): M[i][j]=A[i][j]
        for i in range(7): M[8+i][8+i]=F(b)
        return M
    for X,exp in [(P,diag(I8,0)),(I8,diag(I8,1)),(Z,diag(K,0)),(K,diag(K,1))]:
        if mm(mm(R,diag(X,1)),R)!=exp: bad+=1
    if mm(R,R)!=diag(I8,1) or any(R[i][j]!=R[j][i] for i in range(15) for j in range(15)): bad+=1
print('hub projector identities rho diag(X,I7) rho (P->diag(I,0), H->I, 0->diag(K,0), K->diag(K,I7)) and rho^2=I, rho=rho^T: failures',bad,'of',5*V)
