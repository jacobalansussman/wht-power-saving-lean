import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep import *
import random
def collapse_order(p,keyf):
    res=P(p.k)
    for s in sorted(p.sites(),key=keyf):
        res=res*P(0,{(0,r):v for (si,r),v in p.ops.items() if si==s})
    return res
def PiK(a,keyf): return Cl1(a.m,{g:collapse_order(a.imgs[g],keyf) for g in a.imgs})
rng=random.Random(11)
for name,keyf in [('increasing',lambda s:s),('decreasing',lambda s:-s),('evens-first',lambda s:(s%2,s)),('abs',lambda s:(abs(s),s))]:
    inv=nm=0
    for m in (2,3):
        for it in range(80):
            a=rand_qca(m,rng,rng.randrange(3,25)); b=rand_qca(m,rng,rng.randrange(3,25))
            inv+= (not PiK(a,keyf).valid())
            nm += PiK(a*b,keyf).key()!=(PiK(a,keyf)*PiK(b,keyf)).key()
    print(name,'invalid',inv,'nonmult',nm,'of 160')
