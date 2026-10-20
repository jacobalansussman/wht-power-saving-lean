import sys, os, random; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep import *
def wrap(p,L): # wrap sites mod L (assumes no collisions)
    ops={}
    for (s,r),v in p.ops.items():
        q=(s%L,r); assert q not in ops; ops[q]=v
    return P(p.k,ops)
def psi(a,L):
    imgs={}
    for g,img in a.imgs.items():
        tot=P(0)
        for j in range(L): tot=tot*wrap(img.shift(j),L)
        # must be translation invariant: same letters at every site
        per=[tuple(tot.ops.get((j,r),'I') for r in range(a.m)) for j in range(L)]
        if len(set(per))!=1: return None
        c=sum(1 for v in per[0] if v=='Y')
        sig=P(0,{(0,r):v for r,v in enumerate(per[0]) if v!='I'})
        imgs[g]=P(tot.k+(L-1)*c, sig.ops)
    return Cl1(a.m,imgs)
def conj(c):  # complex conjugate: i -> -i with X, Z fixed: Y -> -Y ; on images in letter form
    imgs={}
    for g,p in c.imgs.items():
        ny=sum(1 for v in p.ops.values() if v=='Y')
        imgs[g]=P((-p.k)+2*ny, p.ops)   # conj(i^k * prod letters) = i^-k * prod conj(letters), conj(Y) = -Y
    return Cl1(c.m,imgs)
rng=random.Random(3); st={'n':0,'L1_eq':0,'L3_eqconj':0,'L1_ne':0,'L3_ne':0,'skip':0}
for m in (1,2,3):
    for it in range(50):
        a=rand_qca(m,rng,rng.randrange(3,20))
        rad=max(max(abs(s) for s in p.sites()) for p in a.imgs.values())
        pi=Pi0(a)
        for L in range(2*rad+3, 2*rad+11, 2):
            ps=psi(a,L)
            if ps is None: st['skip']+=1; continue
            # conjugation also must act on the generators' own phase convention: generators X,Z real -> unchanged
            if L%4==1:
                if ps.key()==pi.key(): st['L1_eq']+=1
                else: st['L1_ne']+=1
            else:
                if ps.key()==conj(pi).key() : st['L3_eqconj']+=1
                else: st['L3_ne']+=1
        st['n']+=1
print(st)
