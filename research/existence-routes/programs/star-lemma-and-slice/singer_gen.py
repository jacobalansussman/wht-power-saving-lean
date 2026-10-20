# Singer (multiplicative-circulant) ansatz, general n, s, prime p (p = 1 mod 2^n - 1). Writes a Singular script.
import sys
n,s,p=map(int,sys.argv[1:4]); N=1<<n; q=N-1
poly={3:0b1011,4:0b10011}[n]
def mul(a,b):
    r=0
    for i in range(n):
        if b>>i&1: r^=a<<i
    for i in range(2*n-2,n-1,-1):
        if r>>i&1: r^=poly<<(i-n)
    return r
pw=[1]
for k in range(q-1): pw.append(mul(pw[-1],2))
assert len(set(pw))==q
def tr(x):
    y=x; t=0
    for _ in range(n): t^=y; y=mul(y,y)
    return t
S0=[k for k in range(q) if tr(pw[k])==0]; S1=[k for k in range(q) if tr(pw[k])==1]
divs=[d for d in range(1,q) if q%d==0]
w=[x for x in range(2,p) if pow(x,q,p)==1 and all(pow(x,d,p)!=1 for d in divs)][0]
iq=pow(q,-1,p)
def cexpr(k):
    return "("+"+".join((f"{iq*pow(w,(-j*k)%q,p)%p}*u({j})" if j else f"{iq*pow(w,(-j*k)%q,p)%p}") for j in range(q))+")"
L=[f"ring r={p},(u(1..{q-1}),v(1..{q-1})),dp;"]
for name,S in (("A",S0),("B",S1)):
    m=len(S); L.append(f"matrix {name}[{m}][{m}];")
    for a in range(m):
        for b in range(m): L.append(f"{name}[{a+1},{b+1}]={cexpr((S[a]-S[b])%q)};")
L.append(f"ideal I=minor(A,{s}),minor(B,{s+1});")   # rank C[S0] <= s-1 (E_00 != 0 takes one), rank C[S1] <= s
L.append(f"for(int j=1;j<={q-1};j++){{I=I,u(j)*v(j)-1;}}")
L.append('ideal G=std(I); print("size,dim:"); print(size(G)); print(dim(G)); if(size(G)<6){print(G);} quit;')
open(f"singer_n{n}_s{s}_p{p}.sing","w").write("\n".join(L)+"\n")
print(n,s,p,"S0",S0,"S1",S1,"root",w)
