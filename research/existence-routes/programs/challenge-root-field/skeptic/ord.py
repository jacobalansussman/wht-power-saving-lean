# independent port of LengthSelection.receipt (skeptic)
import random
def isprime(p):
    if p<2: return False
    i=2
    while i*i<=p:
        if p%i==0: return False
        i+=1
    return True
def nlog2(x): return x.bit_length()-1 if x>0 else 0
def receipt(n):
    stack=nlog2(2*n)+1
    harvest=(16*(stack+4))**2
    oE=nlog2(2*harvest)+2
    if not (2**oE<=n): return n
    # reap (2n) (harvest+1): p = 0..harvest
    vol=1; N=2*n
    p=0
    while p<=harvest:
        if p>2 and vol*p>N: break   # vol fixed from here on: every later p also has vol*p > N
        if p>2 and isprime(p) and vol*p<=N: vol*=p
        p+=1
    q=1
    for _ in range(stack):
        if vol*q<N: q*=2
    return (2*n)*(vol*q)*2**oE
M=2**20*255255
bad=[]
pts=[2**20,2**20+1,2**21-1,2**21,2**26-1,2**26,2**30,2**40-1]+[random.randrange(2**20,2**60) for _ in range(300)]
for n in pts:
    d=receipt(n)
    if d%M or not (0<d<8*n**3) or d%n: bad.append(n)
print('checked',len(pts),'bad',bad)
print('n=2^20-1 ->',receipt(2**20-1)==2**20-1)
# min headroom
print('min 1024n^3/d over pts', min(1024*n**3/receipt(n) for n in pts))
