# Independent skeptic check (own conventions): Pauli strings as per-qubit letters I,X,Y,Z with a global phase in Z4.
# Tests: Lemma 2 (collapse is a valid Clifford map), Lemma 3 (collapse multiplicative), and that gamma = 0,
# on random Clifford QCAs built from elementary layers (one-site H,S,CNOT; row shifts; inter-site CZ; chain maps E).
import random, sys
# single-qubit product table: (a,b) -> (phase exponent mod 4, c), a*b = i^ph c
T = {}
L = 'IXYZ'
for a in L:
    T[('I',a)] = (0,a); T[(a,'I')] = (0,a)
for a in 'XYZ':
    T[(a,a)] = (0,'I')
T[('X','Y')] = (1,'Z'); T[('Y','X')] = (3,'Z')
T[('Y','Z')] = (1,'X'); T[('Z','Y')] = (3,'X')
T[('Z','X')] = (1,'Y'); T[('X','Z')] = (3,'Y')

class P:
    # phase k (i^k), ops: dict (site,row)->letter (non-identity only)
    __slots__=('k','ops')
    def __init__(self,k=0,ops=None):
        self.k=k%4; self.ops=dict(ops or {})
    def __mul__(self,o):
        k=self.k+o.k; ops=dict(self.ops)
        for q,b in o.ops.items():
            a=ops.get(q,'I'); ph,c=T[(a,b)]; k+=ph
            if c=='I': ops.pop(q,None)
            else: ops[q]=c
        return P(k,ops)
    def key(self): return (self.k,tuple(sorted(self.ops.items())))
    def __eq__(self,o): return self.key()==o.key()
    def shift(self,s): return P(self.k,{(si+s,r):v for (si,r),v in self.ops.items()})
    def sites(self): return sorted(set(si for si,_ in self.ops))
    def comm(self,o):  # 0 commute 1 anticommute
        n=0
        for q,a in self.ops.items():
            b=o.ops.get(q)
            if b is not None and b!=a: n+=1
        return n%2
    def sq(self):  # self*self as phase (Paulis letters square to I): i^(2k)
        return (2*self.k)%4

def gen(t,r,site=0): return P(0,{(site,r):t})

class QCA:
    # images of X_r, Z_r at site 0, r in range(m)
    def __init__(self,m,imgs): self.m=m; self.imgs=imgs  # imgs[(t,r)] = P
    def apply(self,p):
        res=P(p.k)
        # p's ops on distinct qubits commute; product of images in any order; Y = i X Z
        for (s,r),a in sorted(p.ops.items()):
            if a=='X': res=res*self.imgs[('X',r)].shift(s)
            elif a=='Z': res=res*self.imgs[('Z',r)].shift(s)
            else:
                res=res*P(1)*self.imgs[('X',r)].shift(s)*self.imgs[('Z',r)].shift(s)
        return res
    def __mul__(self,o):  # (self o other)
        return QCA(self.m,{g:self.apply(o.imgs[g]) for g in o.imgs})

def ident(m): return QCA(m,{(t,r):gen(t,r) for t in 'XZ' for r in range(m)})

def layer_H(m,q):
    d={(t,r):gen(t,r) for t in 'XZ' for r in range(m)}; d[('X',q)]=gen('Z',q); d[('Z',q)]=gen('X',q); return QCA(m,d)
def layer_S(m,q):
    d={(t,r):gen(t,r) for t in 'XZ' for r in range(m)}; d[('X',q)]=gen('Y',q); return QCA(m,d)
def layer_CNOT(m,c,tq):
    d={(t,r):gen(t,r) for t in 'XZ' for r in range(m)}
    d[('X',c)]=gen('X',c)*gen('X',tq); d[('Z',tq)]=gen('Z',c)*gen('Z',tq); return QCA(m,d)
def layer_shift(m,q,s):
    d={(t,r):gen(t,r) for t in 'XZ' for r in range(m)}
    d[('X',q)]=gen('X',q,s); d[('Z',q)]=gen('Z',q,s); return QCA(m,d)
def layer_CZ(m,a,b,s):  # CZ between (a,j) and (b,j+s) for all j, s != 0
    d={(t,r):gen(t,r) for t in 'XZ' for r in range(m)}
    if a==b:
        d[('X',a)]=gen('Z',a,-s)*gen('X',a)*gen('Z',a,s)
    else:
        d[('X',a)]=gen('X',a)*gen('Z',b,s); d[('X',b)]=gen('X',b)*gen('Z',a,-s)
    return QCA(m,d)
def layer_E(m,a,b,s):  # x_a[j] ^= x_b[j-s]: CNOT control (b,j-s) target (a,j); a != b or s != 0
    d={(t,r):gen(t,r) for t in 'XZ' for r in range(m)}
    d[('X',b)]=gen('X',b)*gen('X',a,s); d[('Z',a)]=gen('Z',a)*gen('Z',b,-s); return QCA(m,d)

def collapse(p):
    # ordered product over sites in increasing order of the per-site operators -> one-site Pauli (rows only)
    res=P(p.k)
    for s in p.sites():
        part=P(0,{(0,r):v for (si,r),v in p.ops.items() if si==s})
        res=res*part
    return res

class Cl1:  # one-column Clifford map, images of X_r,Z_r on site 0
    def __init__(self,m,imgs): self.m=m; self.imgs=imgs
    def apply(self,p):
        res=P(p.k)
        for (s,r),a in sorted(p.ops.items()):
            assert s==0
            if a=='X': res=res*self.imgs[('X',r)]
            elif a=='Z': res=res*self.imgs[('Z',r)]
            else: res=res*P(1)*self.imgs[('X',r)]*self.imgs[('Z',r)]
        return res
    def __mul__(self,o): return Cl1(self.m,{g:self.apply(o.imgs[g]) for g in o.imgs})
    def key(self): return tuple((g,self.imgs[g].key()) for g in sorted(self.imgs))
    def valid(self):
        gs=sorted(self.imgs)
        for g in gs:
            if self.imgs[g].sq()!=0: return False
        for i,g in enumerate(gs):
            for h in gs[i+1:]:
                want = 1 if (g[1]==h[1] and g[0]!=h[0]) else 0
                if self.imgs[g].comm(self.imgs[h])!=want: return False
        return True

def Pi0(a): return Cl1(a.m,{g:collapse(a.imgs[g]) for g in a.imgs})

def rand_layer(m,rng):
    c=rng.randrange(6)
    q=rng.randrange(m)
    if c==0: return layer_H(m,q)
    if c==1: return layer_S(m,q)
    if c==2 and m>1:
        t=rng.choice([r for r in range(m) if r!=q]); return layer_CNOT(m,q,t)
    if c==3: return layer_shift(m,q,rng.choice([1,-1]))
    if c==4: return layer_CZ(m,q,rng.randrange(m),rng.choice([1,-1,2]))
    if c==5:
        if m==1: return layer_shift(m,q,1)
        b=rng.choice([r for r in range(m) if r!=q]); s=rng.choice([0,1,-1,2,-2,3])
        return layer_E(m,q,b,s)
    return layer_H(m,q)

def rand_qca(m,rng,nl):
    a=ident(m)
    for _ in range(nl): a=rand_layer(m,rng)*a
    return a

def herm_of(v_ops):
    # Hermitian Pauli from letters (letters I,X,Y,Z products are Hermitian already)
    return P(0,v_ops)

def symbol_parts(a,v):
    # alpha(P_v at site 0) per site: dict site -> per-row letters
    img=a.apply(v)
    return img

def gamma_check(a,m,rng):
    # gamma(v,w) = sum_{d>d'} <M_d v, M_d' w>
    letters='IXYZ'
    v=P(0,{(0,r):rng.choice(letters) for r in range(m)}); v.ops={k:x for k,x in v.ops.items() if x!='I'}
    w=P(0,{(0,r):rng.choice(letters) for r in range(m)}); w.ops={k:x for k,x in w.ops.items() if x!='I'}
    A=a.apply(v); B=a.apply(w)
    sa=A.sites(); sb=B.sites()
    def part(p,s): return P(0,{(0,r):x for (si,r),x in p.ops.items() if si==s})
    g=0
    for d in sa:
        for dp in sb:
            if d>dp: g+=part(A,d).comm(part(B,dp))
    # also direct: collapse homomorphism on v,w: collapse(alpha(v w)) == collapse(alpha v) collapse(alpha w)
    lhs=collapse(a.apply(v*w)); rhs=collapse(A)*collapse(B)
    return g%2, lhs==rhs

def main(seed):
    rng=random.Random(seed)
    stats=dict(qcas=0,pairs=0,invalid=0,nonmult=0,gamma_nonzero=0,hom_fail=0,maxrad=0)
    for m in (1,2,3):
        for it in range(60):
            a=rand_qca(m,rng,rng.randrange(3,25)); b=rand_qca(m,rng,rng.randrange(3,25))
            for x in (a,b):
                stats['qcas']+=1
                if not Pi0(x).valid(): stats['invalid']+=1
                for _ in range(6):
                    g,ok=gamma_check(x,m,rng)
                    stats['gamma_nonzero']+=g; stats['hom_fail']+= (not ok)
                rad=max(max(abs(s) for s in p.sites()) for p in x.imgs.values())
                stats['maxrad']=max(stats['maxrad'],rad)
            stats['pairs']+=1
            if (Pi0(a*b)).key()!=(Pi0(a)*Pi0(b)).key(): stats['nonmult']+=1
    print(seed,stats)
    # sanity controls: collapse must NOT be multiplicative if we use an order that ignores the lag structure? control:
    # reverse-order collapse should also be multiplicative by symmetry; random-order (by hash) collapse should fail sometimes
    return stats

if __name__=='__main__':
    for s in map(int,sys.argv[1:]): main(s)
