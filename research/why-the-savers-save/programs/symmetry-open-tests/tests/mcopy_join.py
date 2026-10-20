#!/usr/bin/env python3
"""t-open test I1 (group 0 of the U8 word, data copied by hand from helpers/d-mcopy-code/group0_spec.out).
Question: the 4 targets' post-cut packet content is 4 Walsh rows of H_8 (rank 4), yet 8 M copies are used.
Information floor if a relay could carry any packet combination = rank of the content restricted to what each
butterfly pair needs.  Frame floor: a relay carrying packets p and q must sit at a frame containing F_p + F_q,
and must meet a target inside that target's end hyperplane.  Which pairs (p,q) have F_p+F_q inside a target's
end frame, and of what rank?"""
import itertools
def basis(vs):
    bb=[]
    for y in vs:
        for z in bb: y=min(y,y^z)
        if y: bb.append(y)
    return bb
def inspan(x,bb):
    for z in bb: x=min(x,x^z)
    return x==0
F={0:[290,162,98,24,6],1:[320,192,24,4,2],2:[320,192,42,24,6],3:[266,138,74,24,6],
   4:[299,171,107,24,6],5:[320,192,19,11,6],6:[320,192,35,24,6],7:[259,131,67,24,6]}
T={169:[320,192,32,16,8,4,2,1],170:[258,130,66,32,16,8,6,1],172:[264,136,72,32,24,4,2,1],175:[258,130,66,32,18,10,6,1]}
W={169:'++++++++',170:'-++--++-',172:'++++----',175:'-+-+-+-+'}
# content rank over Q
import fractions
rows=[[1 if c=='+' else -1 for c in W[t]] for t in T]
def rk(M):
    M=[[fractions.Fraction(x) for x in r] for r in M];r=0
    for c in range(len(M[0])):
        p=next((i for i in range(r,len(M)) if M[i][c]),None)
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]/M[r][c];M[i]=[a-f*b for a,b in zip(M[i],M[r])]
        r+=1
    return r
print('content rank (4 targets x 8 packets):',rk(rows))
for t in T: print('target',t,'end rank',len(basis(T[t])),'packets inside end frame:',[p for p in F if all(inspan(x,basis(T[t])) for x in F[p])])
print('pair joins: dim(F_p+F_q), targets whose end frame contains the join')
cnt=0
for p,q in itertools.combinations(F,2):
    J=basis(F[p]+F[q]); ins=[t for t in T if all(inspan(x,basis(T[t])) for x in J)]
    if ins: cnt+=1
    print(p,q,len(J),ins)
print('pairs whose join fits in some target end frame:',cnt,'of 28')
