# rl-plan: mixed-target moment root (b = 1 - tau, tau root of sum_children n_r r^tau = sum_arrays k_w^tau),
# same formula as rl-vprice/tools/vprice.py and rl-v-logic/partial.py. Controls first, then the variants that
# satisfy Theorem 1G's separation (S): every paid rank < every target rank.
from decimal import Decimal as D, getcontext
from math import log
getcontext().prec = 60
pw = lambda x, t: (D(x).ln() * t).exp()
CALLS = {1: 1560432, 2: 336336, 3: 278724, 4: 84912, 5: 54612, 8: 5328}; W = 66304; HEADS = 11010; TAILS = 11738
def root(children, targets):
    lo, hi = D('0.98'), D(1)
    for _ in range(190):
        x = (lo + hi) / 2
        lhs = sum(n * pw(r, x) for r, n in children.items() if n)
        rhs = sum(n * pw(k, x) for k, n in targets.items() if n)
        if lhs > rhs: lo = x
        else: hi = x
    return 1 - (lo + hi) / 2
def add(d, r, n):
    d = dict(d); d[r] = d.get(r, 0) + n; return d
def sep(children, targets):
    return max(r for r, n in children.items() if n) < min(k for k, n in targets.items() if n)
HUB = D('1.6797549561255445e-3')
def show(name, ch, tg):
    b = root(ch, tg); s = sep(ch, tg)
    rmax = max(r for r, n in ch.items() if n); kmin = min(k for k, n in tg.items() if n)
    th = log(rmax) / log(kmin)
    print(f"{name:58s} b={b:.10e}  vs hub {float((b-HUB)/HUB)*100:+.3f}%  ({float(b-HUB):+.3e})  separated={s} r_max={rmax} k_min={kmin} theta={th:.4f}")
    return b
print("== controls")
show("T74 uniform (expect 1.6797549561e-3)", add(CALLS, 48, HEADS + TAILS), {73: W})
show("T74 all 22,748 truncated (expect 1.8169467780e-3)", CALLS, {73: W - HEADS - TAILS, 25: HEADS + TAILS})
show("T74 heads truncated, tails rank 48 (expect 1.7434773229e-3)", add(CALLS, 48, TAILS), {73: W - HEADS, 25: HEADS})
show("partial f=0.5 t=11374 (expect 1.6983759e-3)", add(CALLS, 24, 2 * (22748 - 11374)), {73: W - 11374, 25: 11374})
print("== separated variants on the computed head count (11,010 heads, 11,738 non-qualifying)")
best = None
for d in range(0, 49):
    tgt = 73 - d                      # head arrays' reduced target
    ch = dict(CALLS)
    if 48 - d: ch = add(ch, 48 - d, HEADS)   # what is left of each head is paid
    cap = tgt - 1                     # largest allowed paid rank
    if 48 <= cap: ch = add(ch, 48, TAILS)
    else:                             # split each tail block into pieces < tgt, as few and as unequal as possible
        rest = 48; pieces = []
        while rest: p = min(cap, rest); pieces.append(p); rest -= p
        for p in pieces: ch = add(ch, p, TAILS)
    if d == 0: tg = {73: W}
    else: tg = {73: W - HEADS, tgt: HEADS}
    b = root(ch, tg)
    ok = sep(ch, tg) if d else None
    if d in (0, 1, 8, 12, 16, 20, 23, 24, 25, 30, 36, 40, 44, 47, 48) or (best and b > best[1]):
        show(f"heads cut by d={d} (target {tgt}), tails {'48' if 48 <= cap else 'split'}", ch, tg)
    if ok and (best is None or b > best[1]): best = (d, b)
print("best separated d:", best[0], f"{best[1]:.10e}")
