#!/usr/bin/env python3
"""Exchange rates for #799's complex side b and PR589's bit side c, computed with pricer.py's engine.
Run: python3 -I levers.py   (a few seconds). Prints a table and writes levers.json next to it.
Exchange rates use the unrounded root (Decimal, 50 digits) of sum_r (r n_r)/(mW) (m/r)^b = 1; any what-if
profile can be certified on the hub's 1e-24 grid with  python3 -I pricer.py --profile ... --m .. --W ..
"""
import sys, json, math, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.environ.get('PRICER_DIR', '.'))  # folder holding pricer.py (it ships with the lever package's tools)
from pricer import CONTROL, C_PR589, kappa_from, price_b, dec
from fractions import Fraction as F
from collections import Counter
from decimal import Decimal as Dd, localcontext

def root(H, m, mW, fallback=0):
    """unrounded b: sum_r r n_r/mW (m/r)^b + fallback m^b = 1 (fallback = 32 m calls/(1e16 W) for the bit surcharge)."""
    with localcontext() as c:
        c.prec = 50
        mW = Dd(mW); terms = [(Dd(r*n)/mW, (Dd(m)/Dd(r)).ln()) for r, n in H.items() if n]
        fb = Dd(F(fallback).numerator)/Dd(F(fallback).denominator); lm = Dd(m).ln()
        lo, hi = Dd(0), Dd('0.03')
        for _ in range(170):
            a = (lo+hi)/2
            if sum(p*(a*l).exp() for p, l in terms) + fb*(a*lm).exp() < 1: lo = a
            else: hi = a
        return (lo+hi)/2
def add(H, d, k=1):
    G = Counter(H)
    for r, n in d.items(): G[r] += k*n
    assert all(v >= 0 for v in G.values()), d
    return {r: n for r, n in G.items() if n}
def mass(H): return sum(r*n for r, n in H.items())
OUT = {}
def row(key, val, sec):
    OUT.setdefault(sec, {})[key] = val; print(f'  {key:78s} {val}')

# ================= complex side (#799) =================
m, W = 73, 7168
H0 = dict(CONTROL['profile']); b0 = root(H0, m, m*W)
LOC = {1:2548, 2:338, 3:430, 4:194, 5:56, 8:8}; IDLE = {1:1, 3:1, 44:1, 48:2396}; CALLS = 72
def B(H, mm=m, mW=m*W): return root(H, mm, mW)
def fmt(x, base=None):
    base = base or b0; return f'{float(x):+.4e} ({float(x/base)*100:+.4f}% )'
S = 'complex'
print('COMPLEX SIDE b (#799): b0 =', f'{float(b0):.13e}')
Fp = sum((r*n/(m*W))*math.log(m/r)*(m/r)**float(b0) for r, n in H0.items())
row('F\'(b0) (mass-weighted log cost incl. growth)', f'{Fp:.6f}', S)
row('mass fraction at rank 48 / at rank<=5', f'{48*2396/(m*W):.4f} / {sum(r*n for r,n in H0.items() if r<=5)/(m*W):.4f}', S)
print(' price list: one extra global child (mass +r, deficit -r)')
for r in (1,2,3,4,5,8,9,44,48):
    row(f'+1 child of rank {r}', fmt(B(add(H0,{r:1}))-b0), S)
print(' deficit / W / m')
row('+1 deficit unit, profile fixed (mW -> mW+1)', fmt(B(H0, m, m*W+1)-b0), S)
row('+1 live family W (=+73 deficit), profile fixed', fmt(B(H0, m, m*(W+1))-b0), S)
row('-1 rank-1 child (=+1 deficit, physical form)', fmt(B(add(H0,{1:-1}))-b0), S)
row('m -> 74 at fixed profile and fixed deficit 2560 (W=(mass+D)/m)', fmt(B(H0, 74, mass(H0)+2560)-b0), S)
row('m -> 72 at fixed profile and fixed deficit 2560', fmt(B(H0, 72, mass(H0)+2560)-b0), S)
row('m -> 74 at fixed W (deficit +7168), profile fixed', fmt(B(H0, 74, 74*W)-b0), S)
print(' helpers (hosting: one unit of helper residual width per call = 72 units moved between rank 48 and the helper run ranks)')
def width_move(r_from_runs, units_per_call=1):
    """remove units_per_call of helper width per call (runs given as {rank: count} per call) and pay
    72*units at rank 48 instead (1.5 rank-48 children per unit): mass-neutral."""
    G = add(H0, r_from_runs, -CALLS)
    extra = CALLS*units_per_call   # mass moved to idle remainder
    q, rem = divmod(extra, 48)
    G = add(G, {48: q}); 
    if rem: G = add(G, {rem: 1})
    return B(G)-b0
for r in (1,2,3,4,5,8):
    row(f'one helper run of rank {r} per call removed, its {r} units hosted->paid at 48', fmt(width_move({r:1}, r)), S)
avg = {1:2548, 2:338, 3:430, 4:194, 5:56}  # non-copy local children
am = mass(avg)
# a "typical" width-w helper: runs in proportion to the local non-copy profile, scaled to mass w (fractional, priced linearly)
def lin_rate(d):  # first-order db for a profile change d (global counts), via F'
    return -sum((r*n/(m*W))*(m/r)**float(b0) for r, n in d.items())/Fp
per_unit_avg = lin_rate({r: CALLS*n/am for r, n in avg.items()}) * -1 + lin_rate({48: CALLS/48})
row('typical helper width unit/call removed (runs ~ local mix), linearised', fmt(Dd(per_unit_avg)), S)
for w in (4,7,8,9):
    row(f'typical helper of width {w} removed (all 72 calls), linearised', fmt(Dd(w*per_unit_avg)), S)
row('helper of width 9 removed, runs all rank 1 (worst) ', fmt(width_move({1:9}, 9)), S)
print(' mass-neutral merges per call (x72), e.g. X4 exchange {2:-4,3:-8,4:+8}')
for d, name in [({1:-2,2:1},'1+1->2'),({1:-1,2:-1,3:1},'1+2->3'),({2:-2,4:1},'2+2->4'),({1:-1,3:-1,4:1},'1+3->4'),({2:-1,3:-1,5:1},'2+3->5'),({4:-2,8:1},'4+4->8'),({2:-4,3:-8,4:8},'X4 joint exchange (#799 vs SW4)')]:
    row(f'merge {name} in every call', fmt(B(add(H0, d, CALLS))-b0), S)
print(' copies (cst = 64 = 8 copy children of rank 8 per call; deficit T[(s-1)v - s cst])')
row('cst 64 -> 63: one copy child 8 -> 7 per call (deficit +72)', fmt(B(add(H0,{8:-1,7:1},CALLS))-b0), S)
row('cst 64 -> 56: one copy child removed per call (deficit +576)', fmt(B(add(H0,{8:-1},CALLS))-b0), S)
print(' rank moves of the idle remainder')
row('one 48 -> 44 (mass -4, deficit +4)', fmt(B(add(H0,{48:-1,44:1}))-b0), S)
row('all 2396 48 -> 44 (deficit +9584)', fmt(B(add(H0,{48:-2396,44:2396}))-b0), S)
row('split one 48 -> 44+4 (mass-neutral, i.e. hosting 4 more units badly)', fmt(B(add(H0,{48:-1,44:1,4:1}))-b0), S)
row('merge the split slot {44,3,1} -> {48} (literal vs cap)', fmt(B(add(H0,{44:-1,3:-1,1:-1,48:1}))-b0), S)
print(' structure what-ifs (same local word, smallest-slot-first hosting model)')
def stages(s, local=LOC, v=112, T=8, a=8, cst=64, h=9):
    K = (s-1)//2; mm = 1+s*a; WW = 2*K*v*T; calls = s*T
    slots = Counter({(2*K-2)*a: 1}) if K>1 else Counter()
    for p in range(1, K):
        for r in (2*p*a, a, (2*K-2*p-2)*a):
            if r: slots[r] += 1
    slots = Counter({r: n*2*v*T for r, n in slots.items()})
    idle_mass = mass(slots); hosted = calls*(mass(local)-2*v*a-cst)
    G = Counter({r: calls*n for r, n in local.items()})
    left = hosted; note = 'all helpers hosted'
    for r in sorted(slots):           # host into the smallest idle slots first (as #799 does: only 48-slots stay paid)
        take = min(slots[r], left//r); slots[r] -= take; left -= take*r
        if 0 < left < r and slots[r]:  # partial slot: the unhosted part is paid as one child
            slots[r] -= 1; G[r-left] += 1; left = 0
        if left == 0: break
    if left > 0: note = f'hosting capacity short by {left} units: they stay unhosted (deficit lowered)'
    for r, n in slots.items():
        if n: G[r] += n
    D = mm*WW - mass(G); assert D == T*((s-1)*v - s*cst) - left, (D,)
    bb = root(dict(G), mm, mm*WW)
    return dict(s=s, m=mm, W=WW, D=D, b=f'{float(bb):.6e}', change=f'{float(bb/b0-1)*100:+.2f}%', maxrank=max(G), note=note)
for s in (7, 9, 11, 13):
    row(f'stages s={s} (m=1+8s), same per-call local word', stages(s), S)
print(' headroom (format-limits-corrected, complex unit, k=1)')
lit = {8: 16128+576+5376, 16: 3584, 32: 3584, 48: 3584}
row('literal R=0 profile at s=9 (no helpers; data climbs as rank-8 runs)', f'{float(root(lit,m,m*W)):.6e}', S)
row('format ceiling X(9) / X(5) / kappa<=', '3.9946e-3 / 5.1863e-3 / 5.1595e-3 (hub, certified)', S)

# ================= bit side (PR589) =================
S = 'bit'
mb, Wb = 64, 28416
HB = {1:381024, 2:400416, 3:22752, 4:26880, 5:29376, 6:17376, 7:26784, 42:288}
BLOC = {1:1323, 2:1342, 3:79, 4:44, 5:102, 6:11, 7:19}; BC = 288
def fbk(H): return F(32*mb*sum(H.values()), 10**16*Wb)
def Bc(H, mW=mb*Wb): return root(H, mb, mW, fbk(H))
c0 = Bc(HB)
print('BIT SIDE c (PR589, with the inherited surcharge): c0 =', f'{float(c0):.13e}', ' certified:', dec(price_b(HB, mb, Wb, True)))
row('surcharge cost (c without surcharge - c)', f'{float(root(HB,mb,mb*Wb)-c0):+.3e}', S)
row('bit idle remainder (global - 288 x local)', str({r: HB[r]-BC*BLOC.get(r,0) for r in HB}), S)
for r in (1,2,7,42):
    row(f'+1 bit child of rank {r}', fmt(Bc(add(HB,{r:1}))-c0, c0), S)
row('+1 bit deficit unit, profile fixed', fmt(Bc(HB, mb*Wb+1)-c0, c0), S)
row('bit cst 63 -> 62 (one copy 7 -> 6 per call, deficit +288)', fmt(Bc(add(HB,{7:-1,6:1},BC))-c0, c0), S)
row('bit: one copy child of rank 7 removed per call (cst 63->56, deficit +2016)', fmt(Bc(add(HB,{7:-1},BC))-c0, c0), S)
# bit helper (width 8 = h, hosted whole): remove its runs per call, the 8x288=2304 freed idle units paid at rank 7 (a) or 42
for runs, name in [({1:8},'runs all rank 1'), ({2:4},'runs rank 2'), ({4:2},'runs rank 4')]:
    G = add(HB, runs, -BC)
    row(f'bit helper removed ({name}); freed 2304 idle units paid at rank 7', fmt(Bc(add(G,{7:2304//7, 2304%7:1}))-c0, c0), S)
    row(f'bit helper removed ({name}); freed 2304 idle units paid at rank 42 (best)', fmt(Bc(add(G,{42:2304//42, 2304%42:1}))-c0, c0), S)
for d, name in [({1:-2,2:1},'1+1->2'),({2:-2,4:1},'2+2->4'),({1:2,3:-1,4:-3,6:1,7:1},'PR589 tune moves (+2e1-e3-3e4+e6+e7)')]:
    row(f'bit merge {name} in every call (x288)', fmt(Bc(add(HB, d, BC))-c0, c0), S)
# kappa cap and targets
kcap = kappa_from(F(2, 1000), C_PR589)[0]
row('kappa cap with PR589 c (b >= c)', dec(kcap, 16), 'kappa')
row('b needed to reach the cap', '>= c = 1.7334270e-3 (+4.147% over b0)', 'kappa')
for target in ('0.00175', '0.0018', '0.002'):
    # smallest min(b,c) giving kappa >= target
    lo, hi = F(17, 10000), F(3, 1000)
    for _ in range(80):
        mid = (lo+hi)/2
        k = kappa_from(mid*(1+F(1,10**12)), mid, check=False)[0]
        if k >= F(target): hi = mid
        else: lo = mid
    row(f'kappa {target} needs min(b,c) >=', f'{float(hi):.7e} (b +{(float(hi)/float(b0)-1)*100:.2f}%, c +{(float(hi)/float(c0)-1)*100:.2f}%)', 'kappa')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'levers.json'), 'w'), indent=1, default=str)
