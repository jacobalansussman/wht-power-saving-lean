from fractions import Fraction as F
from math import comb, log
ok=True
# cap closed form vs verify_hook's proved upper bound c0 = dimN - (n/2)[C(n,k-1) - C(n-2,k-3)]
for m in range(2,9):
    n=2**m
    for k in range(1,n):
        cf=F(k*(n-k)*comb(n,k),2*(n-1))
        c0=k*comb(n,k)-(n//2)*(comb(n,k-1)-(comb(n-2,k-3) if k>=3 else 0))
        if cf!=c0: ok=False; print("MISMATCH",m,k,cf,c0)
print("cap closed form == c0 for m<=8 all k:",ok)
def val(m,k):
    n=2**m; dem=k*comb(n,k)-comb(n,k-1); cap=F(k*(n-k)*comb(n,k),2*(n-1)); return dem,cap,F(dem)/cap
for (m,k) in [(3,2),(3,3),(3,4),(4,3),(4,4),(4,5),(4,6),(4,8),(5,3),(5,4),(5,5),(6,3),(6,4)]:
    d,c,r=val(m,k); n=2**m
    print(m,k,"delta",d,"cap",c,"ratio",r,float(r),"2(n-1)/(n-k+1)=",F(2*(n-1),n-k+1))
for m in [2,3,4,5,6]:
    n=2**m; print("ceiling m",m,F(4*(n-1),n+2),float(F(4*(n-1),n+2)), "smallest cap", log(m/(4*(n-1)/(n+2)))/log(m) if m>2 else 0)
print("m4 k8 actual", F(86592,27456), "cap(4,8)", val(4,8)[1])
for m,L in [(3,F(14,5)),(4,F(41,13)),(5,F(31,14)),(6,F(63,31))]:
    print("moves-only cap m",m,L,log(m/L)/log(m))
