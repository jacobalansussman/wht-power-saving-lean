"""v-s-fix/algebra_check.py: exact check of s-fix's escape-(b) algebra, in s-fix's model and in a fuller model
where a target may also read OTHER targets (legal in the local reading: every P-bank target ends at local H).
Variables: x_old, y_old independent; z = x_old + y_old.  Source T at stage i+1 holds c_T = z_T + (E_z z + E_y y_old)_T.
Target S starts with x_old_S and must end with -y_old_S, i.e. receive exactly -z_S.
  s-fix model:  receive = M c                          -> M(I+E_z) = -I, M E_y = 0   (sign convention -I here)
  fuller model: receive = M c + N x_old, N[S,S] = 0     -> M(I+E_z+E_y) = -I, N = M E_y, diag(N) = 0
usage: python3 -I algebra_check.py"""
from fractions import Fraction as F
import random
random.seed(1)
def mat(n,f): return [[F(f(i,j)) for j in range(n)] for i in range(n)]
def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def add(A,B,s=1): return [[a+s*b for a,b in zip(r,q)] for r,q in zip(A,B)]
def inv(A):
    n=len(A); M=[r[:]+[F(int(i==j)) for j in range(n)] for i,r in enumerate(A)]
    for c in range(n):
        p=next(r for r in range(c,n) if M[r][c]!=0); M[c],M[p]=M[p],M[c]
        iv=1/M[c][c]; M[c]=[x*iv for x in M[c]]
        for r in range(n):
            if r!=c and M[r][c]!=0: f=M[r][c]; M[r]=[x-f*y for x,y in zip(M[r],M[c])]
    return [r[n:] for r in M]
def mv(A,v): return [sum(a*b for a,b in zip(r,v)) for r in A]
def simulate(M,N,Ez,Ey,trials=20):
    n=len(M); ok=True
    for _ in range(trials):
        xo=[F(random.randint(-9,9)) for _ in range(n)]; yo=[F(random.randint(-9,9)) for _ in range(n)]
        z=[a+b for a,b in zip(xo,yo)]
        c=[a+b+d for a,b,d in zip(z,mv(Ez,z),mv(Ey,yo))]
        final=[xo[s]+r1+r2 for s,(r1,r2) in enumerate(zip(mv(M,c),mv(N,xo)))]
        ok&= final==[-y for y in yo]
    return ok
I2=mat(2,lambda i,j:i==j); Z2=mat(2,lambda i,j:0)
# (1) s-fix's own model: M = -(I+E_z)^{-1}; any E_y != 0 contradicts M E_y = 0 since M is invertible (any number of carriers)
Ez=mat(3,lambda i,j:0); Ez[0]=[F(0),F(2),F(-1)]          # single carrier T0=0, s = (0,2,-1)
M=[[ -x for x in r] for r in inv(add(mat(3,lambda i,j:i==j),Ez))]
print('(1) s-fix model, single carrier E_z=e0 s^T, s=(0,2,-1): diag(M) =',[str(M[i][i]) for i in range(3)],
      '; simulate (N=0, E_y=0):',simulate(M,mat(3,lambda i,j:0),Ez,mat(3,lambda i,j:0)))
# (2) fuller model: carrier 0 holds r.y_old with r = e_1 (pure y_old content, E_z = 0)
Ey=[[F(0),F(1)],[F(0),F(0)]]
M=[[-x for x in r] for r in inv(add(I2,Ey))]; N=mm(M,Ey)
print('(2) fuller model, carrier 0 holds y_old_1: M =',[[str(x) for x in r] for r in M],' N =',[[str(x) for x in r] for r in N],
      ' diag N zero:',all(N[i][i]==0 for i in range(2)),' simulate:',simulate(M,N,Z2,Ey))
print('    -> target 0 reads source 1 AND the carried y_old_1 in source 0, and cancels x_old_1 from target 1:'
      ' the carried y_old IS read; s-fix\'s "M E_y = 0, can never be read" holds only without target->target reads.')
print('    diag(M) =',[str(M[i][i]) for i in range(2)],'-> the diagonal is still nonzero (no saving).')
# (3) fuller model, single carrier, random w: diag(M) never zero (rank-one update), N = M E_y with diag N = 0 when admissible
bad=0; tried=0; adm=0
for _ in range(2000):
    n=random.randint(2,5); T0=random.randrange(n)
    s=[F(random.randint(-2,2)) for _ in range(n)]; r=[F(random.randint(-2,2)) for _ in range(n)]
    Ez=mat(n,lambda i,j:0); Ey=mat(n,lambda i,j:0); Ez[T0]=s[:]; Ey[T0]=r[:]
    A=add(add(mat(n,lambda i,j:i==j),Ez),Ey)
    if 1+s[T0]+r[T0]==0: continue
    tried+=1; M=[[-x for x in row] for row in inv(A)]; N=mm(M,Ey)
    if any(N[i][i]!=0 for i in range(n)): continue
    adm+=1
    if any(M[i][i]==0 for i in range(n)): bad+=1
    assert simulate(M,N,Ez,Ey,3)
print('(3) fuller model, single carrier (random s, r, n<=5): %d invertible, %d admissible (diag N = 0), zero diagonal in M: %d' % (tried,adm,bad))
# (4) s-fix's permutation door, both models: E_z swaps z_0, z_1 (E_y = 0)
Ez=[[F(-1),F(1)],[F(1),F(-1)]]
M=[[-x for x in r] for r in inv(add(I2,Ez))]
print('(4) permutation door E_z = swap: diag(M) =',[str(M[i][i]) for i in range(2)],' simulate:',simulate(M,Z2,Ez,Z2))
