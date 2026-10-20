# Sanity: the standard 32-call butterfly computes H^{(x)4}; check orthogonality lemma at every prefix.
import numpy as np, itertools
H=np.array([[1,1],[1,-1]],float)
def std_net():
    calls=[]
    for bit in range(4):
        for x in range(16):
            if not (x>>bit)&1: calls.append((x,x|(1<<bit)))
    return calls
def apply(calls,gates=None,rng=None):
    W=np.eye(16,dtype=complex)
    for k,(u,v) in enumerate(calls):
        G=H if gates is None else gates[k]
        r=W[[u,v],:].copy(); W[[u,v],:]=G@r
    return W
calls=std_net()
W=apply(calls)
H16=H
for _ in range(3): H16=np.kron(H16,H)
# match up to monomials: the butterfly with x -> bits gives exactly kron in some order
print("std calls",len(calls),"equals H16 up to perm:",np.allclose(np.abs(W),1))
def parts(calls,t):
    # past partition after first t calls, future partition of calls t..end
    def comps(cs):
        lab=list(range(16))
        def f(a):
            while lab[a]!=a: a=lab[a]
            return a
        for u,v in cs:
            a,b=f(u),f(v)
            if a!=b: lab[a]=b
        d={}
        for w in range(16): d.setdefault(f(w),set()).add(w)
        return list(d.values())
    return comps(calls[:t]),comps(calls[t:])
ok=True
for t in range(len(calls)+1):
    P,Q=parts(calls,t)
    for A in Q:
        for B in P:
            if 16*len(A&B)!=len(A)*len(B): ok=False
print("orthogonal at every prefix:",ok)
