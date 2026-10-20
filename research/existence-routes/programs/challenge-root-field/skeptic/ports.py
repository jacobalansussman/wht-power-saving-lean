def nlog2(x): return x.bit_length()-1 if x>0 else 0
def ports(n):
    stack=nlog2(2*n)+1; h=(16*(stack+4))**2; oE=nlog2(2*h)+2
    return 2**oE<=n
bad=[n for n in range(1,2**20) if ports(n)]
print('Ports true below 2^20:',len(bad), bad[:5])
# monotone above: check every power-of-two boundary up to 2^200 where stackN changes
import itertools
viol=[k for k in range(20,200) if not (ports(2**k) and ports(2**(k+1)-1))]
print('Ports fails at some boundary 2^k..2^(k+1)-1 for k in [20,200):',viol)
