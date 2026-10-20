# Skeptic check (exact, pure Python, < 1 s, < 20 MB).
# (a) Bluestein analogue for Z_2^n: (-1)^{|x & y|} = i^{|x|} i^{|y|} i^{-|x ^ y|}, so H_n = D K D with
#     K[x][y] = (-i)^{|x^y|} an XOR-convolution (dyadic convolution) by a fixed kernel and D = diag(i^{|x|}).
#     Hence L(WHT) <= 3 C_xorconv + 2N (linearisation) and C_xorconv <= 3 L(WHT) + N: same log-exponent.
# (b) W^{(x)n} (W = e100+e010+e001) equals the subset-convolution tensor sum_{A cap B = 0} e_A (x) e_B (x) e_{A cup B}
#     after relabelling the third index by complement, checked as supports for n = 1..4.
from itertools import product
pc=lambda v: bin(v).count('1')
I=[1,1j,-1,-1j]
for n in range(1,7):
    N=1<<n
    ok=all(((-1)**pc(x&y)) == I[pc(x)%4]*I[pc(y)%4]*I[(-pc(x^y))%4] for x in range(N) for y in range(N))
    print("n=%d: (-1)^{|x&y|} == i^|x| i^|y| i^-|x^y| for all x,y: %s"%(n,ok))
for n in range(1,5):
    N=1<<n; full=N-1
    W=set()
    for trip in product([(1,0,0),(0,1,0),(0,0,1)],repeat=n):
        a=sum(t[0]<<k for k,t in enumerate(trip)); b=sum(t[1]<<k for k,t in enumerate(trip)); c=sum(t[2]<<k for k,t in enumerate(trip))
        W.add((a,b,c))
    SC={(A,B,full^(A|B)) for A in range(N) for B in range(N) if A&B==0}
    print("n=%d: support of W^(x)n == subset-convolution support (third index complemented): %s (%d terms)"%(n,W==SC,len(W)))
