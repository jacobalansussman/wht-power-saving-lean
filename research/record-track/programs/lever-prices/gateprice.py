# hm-gate: own mixed-target moment root, written from the formula (not copied from rl-plan/rescue.py).
# tau solves sum_children n_r r^tau = sum_arrays k_w^tau ; b = 1 - tau. Uniform case = hub moment equation
# sum_r (r n_r)/(m W) (m/r)^b = 1 (multiply through by W m^(1-b)).
import sys, json
from decimal import Decimal as D, getcontext
getcontext().prec = 80
LN = {}
def ln(x):
    if x not in LN: LN[x] = D(x).ln()
    return LN[x]
def f(children, targets, t):
    return sum(n * (ln(r) * t).exp() for r, n in children.items() if n) - sum(n * (ln(k) * t).exp() for k, n in targets.items() if n)
def root(children, targets):
    lo, hi = D('0.9'), D(1)          # f(lo) < 0 < f(hi) when there is a saving
    assert f(children, targets, hi) < 0 or True
    for _ in range(260):
        x = (lo + hi) / 2
        if f(children, targets, x) > 0: lo = x   # children dominate at small t? check sign below
        else: hi = x
    return 1 - (lo + hi) / 2
CALLS = {1: 1560432, 2: 336336, 3: 278724, 4: 84912, 5: 54612, 8: 5328}
W = 66304; IDLE = 22748; HUB = D('0.001679754956125544521335')
def merge(*ds):
    o = {}
    for d in ds:
        for k, v in d.items(): o[k] = o.get(k, 0) + v
    return o
def sep(ch, tg):
    return max(r for r, n in ch.items() if n) < min(k for k, n in tg.items() if n)
out = {}
def show(name, ch, tg):
    # sign convention check: at t=1 children mass < target mass (deficit > 0) => f(1) < 0; at t=0 counts children > arrays => f(0) > 0
    assert f(ch, tg, D(1)) < 0 and f(ch, tg, D(0)) > 0, name
    b = root(ch, tg)
    rel = (b - HUB) / HUB * 100
    out[name] = {"b": f"{b:.13e}", "vs_hub_percent": f"{rel:+.4f}", "separated": sep(ch, tg)}
    print(f"{name:66s} b={b:.10e} {rel:+.3f}% sep={sep(ch, tg)}")
    return b
show("control uniform T74", merge(CALLS, {48: IDLE}), {73: W})
for H in (11010, 9443, 12000, 13000, 14000, 15000, 16576):
    T = IDLE - H
    show(f"full cut, heads={H}, others 24+24", merge(CALLS, {24: 2 * T}), {73: W - H, 25: H})
    show(f"half cut, heads={H} (keep 24, target 49), others 48", merge(CALLS, {48: T, 24: H}), {73: W - H, 49: H})
json.dump(out, open(sys.argv[1], 'w'), indent=1)
