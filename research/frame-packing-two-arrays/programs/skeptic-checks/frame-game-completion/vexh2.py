"""vexh2.py (v-pk2-complete): EXHAUSTIVE homomorphism check at block 2, exact in C[G].
BFS over Sp(4,2) from legitimate letters gives one lift g_lambda for each of the 720 elements. For every lambda and
every generator x (all dir z, phase letters, elementary address maps, a shift times i):
   Lambda(xbar) Lambda(lambda) == Lambda(xbar lambda)   with Lambda(xbar lambda) computed from the STORED lift of
   xbar lambda (a different matrix from x g_lambda in general: this also checks independence of the lift).
Since the generators generate Sp and Lambda(1) = 1, this proves (exactly, by computation) that Lambda is a unital
homomorphism at block 2."""
import numpy as np, time
from veng import *
E = VEng(2); n = E.n; t0 = time.time()
gens = [E.dirm(z) for z in range(1, n)] + [E.phase_letter(z, c) for z in range(n) for c in (0, 1)] + \
       [E.addr_mat([1, 3]), E.addr_mat([3, 2]), E.addr_mat([2, 1]), 1j * E._X(3)]
lift = {E.skey(E.gbar(E.I)): E.I}; fr = [E.I]
while fr:
    g = fr.pop()
    for x in gens:
        h = x @ g; k = E.skey(E.gbar(h))
        if k not in lift: lift[k] = h; fr.append(h)
print("Sp(4,2) elements reached:", len(lift), flush=True)
Dc = {}
def Lam(g):
    k = E.key(g)
    if k not in Dc:
        out = {}
        for (u, j), c in E.D(g).items():
            kk = E.key(g @ E.nmat(u, j)); out[kk] = out.get(kk, GQ(0)) + c
        Dc[k] = {a: b for a, b in out.items() if not b.iszero()}
    return Dc[k]
def cmul(X, Y):
    out = {}
    for kx, vx in X.items():
        Mx = E.store[kx]
        for ky, vy in Y.items():
            k = E.key(Mx @ E.store[ky]); out[k] = out.get(k, GQ(0)) + vx * vy
    return {k: (v.r, v.i) for k, v in out.items() if not v.iszero()}
ok = bad = 0
for k, g in lift.items():
    Lg = Lam(g)
    for x in gens:
        lhs = cmul(Lam(x), Lg)
        rhs = {a: (b.r, b.i) for a, b in Lam(lift[E.skey(E.gbar(x @ g))]).items()}
        if lhs == rhs: ok += 1
        else: bad += 1
print(f"block 2 exhaustive: {ok} pairs (generator, element) satisfy Lambda(x)Lambda(l) = Lambda(x l) exactly, {bad} fail; "
      f"Lambda(1) = 1: {Lam(E.I) == {E.key(E.I): GQ(1)} or {a:(b.r,b.i) for a,b in Lam(E.I).items()} == {E.key(E.I):(1,0)}}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
