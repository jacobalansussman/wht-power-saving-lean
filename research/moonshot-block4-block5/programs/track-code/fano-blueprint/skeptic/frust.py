"""v-ms-fano: independent check of Lemma P(d) and the RM(1,3) coset claims.
Frustration tested by linear algebra (not by enumerating characters): pieces x1 + lam_j Z^{v_j} x2 separate iff
there are z_j in Q[F_2^3] with sum z_j = 1 and sum lam_j z_j Z^{v_j} = 0. Solve this exact linear system over Q
(8 unknowns per piece) with sympy; compare with the affine-function test."""
import itertools
from fractions import Fraction as F
def rank(M):
    M=[[F(x) for x in r] for r in M]; r=0; nc=len(M[0]) if M else 0
    for c in range(nc):
        p=next((i for i in range(r,len(M)) if M[i][c]!=0),None)
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        for i in range(len(M)):
            if i!=r and M[i][c]!=0:
                f=M[i][c]/M[r][c]; M[i]=[a-f*b for a,b in zip(M[i],M[r])]
        r+=1
    return r

pts = list(itertools.product((0,1),repeat=3))
idx = {p:i for i,p in enumerate(pts)}
def add(a,b): return tuple((x+y)%2 for x,y in zip(a,b))
def sep_linear(pattern):
    n=len(pattern); A=[]; b=[]
    # z_j = sum_h z_{j,h} h ; coefficient of g in sum_j z_j is sum_j z_{j,g}; in sum_j lam_j z_j Z^{v_j} it is sum_j lam_j z_{j,g+v_j}
    for g in pts:
        r1=[0]*(8*n); r2=[0]*(8*n)
        for j in range(n):
            r1[8*j+idx[g]]+=1
            r2[8*j+idx[add(g,pattern[j][0])]]+=pattern[j][1]
        A+=[r1,r2]; b+=[1 if g==(0,0,0) else 0, 0]
    aug=[r+[bb] for r,bb in zip(A,b)]
    return rank(A)==rank(aug)
def affine(pattern):
    for a in pts:
        vals={l*(-1)**(sum(x*y for x,y in zip(a,v))%2) for v,l in pattern}
        if len(vals)==1: return True
    return False
bad=0; tot=0; minsize=None
# every subset S of F_2^3 and every sign pattern: linear separability == not affine
for r in range(1,9):
    for S in itertools.combinations(pts,r):
        for signs in itertools.product((1,-1),repeat=r):
            if signs[0]==-1: continue  # global sign irrelevant
            pat=list(zip(S,signs)); tot+=1
            ls=sep_linear(pat); af=affine(pat)
            if ls==af: bad+=1
            if ls and minsize is None: minsize=r; print('smallest frustrated (distinct points):',r,pat)
    print('size',r,'checked so far',tot,'mismatches',bad,flush=True)
# lines never frustrated
lines=[(u,v,add(u,v)) for u,v in itertools.combinations(pts[1:],2)]
print('lines frustrated ever:', any(sep_linear(list(zip(L,s))) for L in lines for s in itertools.product((1,-1),repeat=3)))
# affine planes: frustrated iff odd # minus
ok=True
for u,v in itertools.combinations(pts[1:],2):
    for base in pts:
        P=[base,add(base,u),add(base,v),add(base,add(u,v))]
        if len(set(P))<4: continue
        for s in itertools.product((1,-1),repeat=4):
            if sep_linear(list(zip(P,s))) != (s.count(-1)%2==1): ok=False
print('affine planes: frustrated iff odd number of minus:',ok)
# cosets of RM(1,3)
code={tuple((c+sum(x*y for x,y in zip(a,v)))%2 for v in pts) for a in pts for c in (0,1)}
cos=set(); leaders=[]
for w in itertools.product((0,1),repeat=8):
    C=frozenset(tuple((x+y)%2 for x,y in zip(w,c)) for c in code)
    if C in cos: continue
    cos.add(C); leaders.append(min(sum(x) for x in C))
print('cosets',len(cos),'leader weights',sorted(leaders),'covering radius',max(leaders))
# weight-k code words: check [8,4,4]
from collections import Counter
print('code weight distribution',sorted(Counter(sum(c) for c in code).items()))
