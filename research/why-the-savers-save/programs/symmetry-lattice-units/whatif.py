"""What-if prices in the hub's format (pricer.py = exact hub moment engine), calibrated on the #1042 word.

Model of one word (k = 1, s stages, T replicas, v ports, window h, a = h - 1, m = 1 + s a):
  W = T v (s - 1),  D = T[(s - 1) v - s cst],  rank mass = m W - D  (the hub's (Dfc)/(Win) identities; #1042: 73*14336-5120).
  per call (T s calls): data climbs = 2v registers climbing a each, cut into runs with #1042's per-port run shape
     (x and y classes from k8_census, mass 8 per port each) scaled to mass a; copies = the decoder's copy children;
     helper runs = #1042's per-port helper-run shape (3730 width units per 112 ports), times a helper factor.
  idle: the B_K schedule's idle slots (2 per port per replica, ranks idle_slots(K, a)), helper mass hosted in them
     smallest slot first (hub (Hst)), what is left paid as one child per slot (the hub's literal remainder).
Calibration: with #1042's own unit (v 112, h 9, cst 64) this model gives the b printed as 'model #1042' (hub 1.674989e-3).
"""
import sys, os, json
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hubcopy'))
from pricer import price_b, kappa_from, C_PR589
from skel import idle_slots

X_RUNS = {1: 406, 2: 44, 3: 30, 4: 58, 5: 16}      # #1042 (k8 word) per call, 112 ports, mass 896
Y_RUNS = {1: 331, 2: 102, 3: 37, 4: 30, 5: 26}     # mass 896
H_RUNS = {1: 1675, 2: 272, 3: 365, 4: 54, 5: 40}   # helper runs, mass 3730, 480 helpers

def scaled(shape, factor):
    """run shape (rank -> count) with every count multiplied by factor (a Fraction); returned as Fractions."""
    return {r: F(n) * factor for r, n in shape.items()}

def word(v, h, copy_children, s=9, T=16, helper_factor=F(1), data_shape=True):
    a = h - 1; m = 1 + s * a; K = (s - 1) // 2
    cst = sum(copy_children)
    W = T * v * (s - 1); D = T * ((s - 1) * v - s * cst)
    calls = T * s
    prof = {}
    def add(r, n):
        if n: prof[r] = prof.get(r, 0) + n
    # data runs
    for shape in (X_RUNS, Y_RUNS):
        if data_shape:
            sc = scaled(shape, F(v, 112) * F(a, 8) * calls)
            for r, n in sc.items(): add(r, n)
        else:
            add(a, 2 * 0 + v * calls)              # literal: one child of rank a per register per stage
    for r in copy_children: add(r, calls)
    hs = {r: int(round(n)) for r, n in scaled(H_RUNS, F(v, 112) * helper_factor * calls).items()}
    for r, n in hs.items(): add(r, n)
    hmass = sum(r * n for r, n in hs.items())
    # idle slots, hosting smallest first (over the whole multiset of slots of all ports and replicas)
    slots = sorted(idle_slots(K, a) * (2 * v * T))
    left = hmass
    for r in slots:
        if left >= r: left -= r; continue
        if left > 0: add(r - left, 1); left = 0; continue
        add(r, 1)
    assert left == 0, 'helpers exceed the idle capacity'
    # integerise: round counts, then put the mass error on rank 1 (tiny)
    mass_target = m * W - D
    H = {r: int(round(n)) for r, n in prof.items()}
    err = mass_target - sum(r * n for r, n in H.items())
    H[1] = H.get(1, 0) + err
    assert H[1] >= 0
    H = {int(r): n for r, n in H.items() if n and r < m}
    # children of rank >= m cost nothing in the engine; drop them (their mass is still counted in m W)
    return H, m, W, D, cst

def price(name, v, h, copies, side, **kw):
    H, m, W, D, cst = word(v, h, copies, **kw)
    # pricer asserts m^2 > 2 max(r)^2; children above m/sqrt2 would need halving: report instead
    mx = max(H)
    ok = m * m > 2 * mx * mx
    if not ok:
        return dict(name=name, error='halving constraint fails (max child %d, m %d)' % (mx, m))
    b = price_b(H, m, W, bad=(side == 'bit'), certify=True)
    rec = dict(name=name, side=side, v=v, h=h, cst=cst, x=float(F(cst, v)), s=kw.get('s', 9), T=kw.get('T', 16), m=m, W=W, D=D,
               helper_factor=str(kw.get('helper_factor', 1)), helpers_per_call=round(480 * v / 112 * float(kw.get('helper_factor', 1))),
               saving=float(b), saving_str='%.7e' % float(b))
    return rec

if __name__ == '__main__':
    rows = []
    rows.append(price('model #1042 (E8 P112, calibration)', 112, 9, [8] * 8, 'complex'))
    rows.append(price('E8 P111 bit in the same model', 111, 8, [7] * 9, 'bit'))
    for hf in (F(1), F(15, 8)):
        rows.append(price('BW16 bit, ovoid decoder, helpers x%s per port' % hf, 2160, 16, [8] * 162, 'bit', helper_factor=hf))
    rows.append(price('BW16 bit, ovoid decoder, NO helpers', 2160, 16, [8] * 162, 'bit', helper_factor=F(0)))
    rows.append(price('BW16 bit, hypothetical x=0.533 (72 totals of rank 16), helpers x1', 2160, 16, [16] * 72, 'bit'))
    rows.append(price('BW16 complex, h=72, cst=0, helpers x1 (best imaginable)', 2160, 72, [], 'complex'))
    for r in rows: print(json.dumps(r))
    json.dump(rows, open('whatif.json', 'w'), indent=1)
