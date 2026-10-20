# Character blocks of the block-3 hook demand (m=3, f=5), exact over Z[i], and mod-p certificates.
# Row reps tau = (0, u1<u2<u3<u4), u_i nonzero; col reps sigma = (0, 0, t1<t2<t3), t_i nonzero.
# Ghat_s[r,c] = sum_u chi_s(u) * Pr(u, 0) * det_{b,a}( Pr(u_b ^ u, t_a) ),  t_0 = 0, (t_1,t_2,t_3) = rho
# Pr(x, y) = (-i)^{popcount(x^y)}.  (Gaussian integers as Python complex with int parts -> use tuples)
import itertools, sys, json
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gadd(a,b): return (a[0]+b[0], a[1]+b[1])
def gneg(a): return (-a[0], -a[1])
MI = [(1,0),(0,-1),(-1,0),(0,1)]  # (-i)^k
def pc(x): return bin(x).count('1')
def Pr(x,y): return MI[pc(x^y)%4]
def det4(Mx):
    # Laplace on rows 0,1 via 2x2 minors (matches the Lean det4)
    def m(r0,r1,i,j): return gadd(gmul(Mx[r0][i],Mx[r1][j]), gneg(gmul(Mx[r0][j],Mx[r1][i])))
    tot=(0,0)
    for (i,j) in itertools.combinations(range(4),2):
        k,l=[c for c in range(4) if c not in (i,j)]
        sgn = (-1)**(i+j+1)  # sign of the 2x2 Laplace term: (-1)^{(0+1)+(i+j)} with rows 0,1 -> (-1)^(1+i+j)
        term = gmul(m(0,1,i,j), m(2,3,k,l))
        if sgn<0: term=gneg(term)
        tot=gadd(tot,term)
    return tot
def det_perm(Mx):
    tot=(0,0)
    for p in itertools.permutations(range(4)):
        s=1
        for i in range(4):
            for j in range(i+1,4):
                if p[i]>p[j]: s=-s
        t=(1,0)
        for i in range(4): t=gmul(t,Mx[p[i]][i])
        tot=gadd(tot, t if s>0 else gneg(t))
    return tot
rows=[ (0,)+U for U in itertools.combinations(range(1,8),4)]
cols=[ (0,0)+T for T in itertools.combinations(range(1,8),3)]
def dot(s,u): return pc(s&u)%2
def entry(s,tau,sig):
    tot=(0,0)
    for u in range(8):
        Mx=[[Pr(tau[b+1]^u, sig[a+1]) for a in range(4)] for b in range(4)]
        d=det4(Mx)
        assert d==det_perm(Mx)
        t=gmul(Pr(tau[0]^u, sig[0]), d)
        if dot(s,u): t=gneg(t)
        tot=gadd(tot,t)
    return tot
p=int(sys.argv[1]) if len(sys.argv)>1 else 65537
I=[x for x in range(2,p) if pow(x,2,p)==p-1][0] if p<10**6 else None
if I is None:
    # find sqrt(-1) mod p for p = 1 mod 4
    for a in range(2,1000):
        c=pow(a,(p-1)//4,p)
        if c*c%p==p-1: I=c; break
print("p",p,"I",I, file=sys.stderr)
def modp(g): return (g[0]+g[1]*I)%p
def rank_select(A):
    # greedy row/col selection via elimination mod p; returns pivot rows, pivot cols
    A=[r[:] for r in A]; n=len(A); m=len(A[0])
    piv_r=[]; piv_c=[]
    R=[r[:] for r in A]
    used=[False]*n
    for c in range(m):
        pr=None
        for r in range(n):
            if not used[r] and R[r][c]%p: pr=r; break
        if pr is None: continue
        used[pr]=True; piv_r.append(pr); piv_c.append(c)
        inv=pow(R[pr][c],p-2,p)
        for r in range(n):
            if r!=pr and R[r][c]%p:
                f=R[r][c]*inv%p
                R[r]=[(x-f*y)%p for x,y in zip(R[r],R[pr])]
    return piv_r,piv_c
def inv_mat(A):
    n=len(A); M=[A[i][:]+[1 if i==j else 0 for j in range(n)] for i in range(n)]
    for c in range(n):
        pr=next(r for r in range(c,n) if M[r][c]%p)
        M[c],M[pr]=M[pr],M[c]
        inv=pow(M[c][c],p-2,p); M[c]=[x*inv%p for x in M[c]]
        for r in range(n):
            if r!=c and M[r][c]:
                f=M[r][c]; M[r]=[(x-f*y)%p for x,y in zip(M[r],M[c])]
    return [row[n:] for row in M]
out={"p":p,"I":I,"rows":rows,"cols":cols,"blocks":[]}
total=0
for s in range(8):
    G=[[entry(s,t,c) for c in cols] for t in rows]
    Gp=[[modp(g) for g in row] for row in G]
    pr,pcs=rank_select(Gp)
    total+=len(pr)
    sub=[[Gp[r][c] for c in pcs] for r in pr]
    B=inv_mat(sub) if len(pr)==len(pcs) else None
    # verify
    n=len(pr)
    ok=all(sum(B[i][k]*sub[k][j] for k in range(n))%p==(1 if i==j else 0) for i in range(n) for j in range(n))
    print("s",s,"rank",len(pr),"inverse ok",ok, "max|G|", max(abs(g[0])+abs(g[1]) for row in G for g in row), file=sys.stderr)
    out["blocks"].append({"s":s,"rsel":pr,"csel":pcs,"B":B,"G":[[G[r][c] for c in pcs] for r in pr]})
print("total",total,file=sys.stderr)
json.dump(out,open(sys.argv[2] if len(sys.argv)>2 else "cert.json","w"))
