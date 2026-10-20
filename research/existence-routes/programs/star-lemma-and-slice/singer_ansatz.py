# Singer (multiplicative-circulant) ansatz for the local AG(4,2), s=3 problem; writes a Singular script over F_31.
# Points = F_16 (alpha^4 = alpha + 1). E commutes with x -> alpha x: E = E_00 (+) C, C circulant on Z_15 (k = log x).
# Hyperplane orbits: H0 = {Tr x = 0} (contains 0) and H1 = {Tr x = 1}. Need rank C[S0,S0] <= 2 and rank C[S1,S1] <= 3,
# C nonsingular (all 15 eigenvalues nonzero).  Unknowns: eigenvalues u_j (u_0 = 1 by scaling) and inverses v_j.
import itertools
def mul(a,b):
    r=0
    for i in range(4):
        if b>>i&1: r^=a<<i
    for i in range(7,3,-1):
        if r>>i&1: r^=0b10011<<(i-4)
    return r
pw=[1]
for k in range(14): pw.append(mul(pw[-1],2))
assert len(set(pw))==15
def tr(x):
    y=x; s=0
    for _ in range(4): s^=y; y=mul(y,y)
    assert s in (0,1); return s
S0=[k for k in range(15) if tr(pw[k])==0]; S1=[k for k in range(15) if tr(pw[k])==1]
print("S0",S0,"S1",S1)
p=31
w=[x for x in range(2,p) if pow(x,15,p)==1 and all(pow(x,d,p)!=1 for d in (1,3,5))][0]
inv15=pow(15,-1,p)
# c_k = (1/15) sum_j u_j w^{-jk}
def cexpr(k):
    terms=[]
    for j in range(15):
        coef=inv15*pow(w,(-j*k)%15,p)%p
        terms.append(f"{coef}*u({j})" if j>0 else f"{coef}")
    return "("+"+".join(terms)+")"
lines=[f"ring r=31,(u(1..14),v(1..14)),dp;","option(redSB);"]
# matrices
def mat(name,S):
    n=len(S)
    lines.append(f"matrix {name}[{n}][{n}];")
    for a in range(n):
        for b in range(n):
            lines.append(f"{name}[{a+1},{b+1}]={cexpr((S[a]-S[b])%15)};")
mat("A",S0); mat("B",S1)
lines.append("ideal I=minor(A,3),minor(B,4);")
lines.append("for(int j=1;j<=14;j++){I=I,u(j)*v(j)-1;}")
lines.append("int t=timer; ideal G=std(I); print(\"dim/size:\"); print(size(G)); print(G[1]); print(dim(G)); print(timer-t);")
lines.append("quit;")
open("singer_ansatz.sing","w").write("\n".join(lines)+"\n")
print("primitive 15th root mod 31:",w)
