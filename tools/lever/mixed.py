# lever-reprice: mixed-target moment root, the same formula as rl-plan/tools/rescue.py (b = 1 - tau, tau the root of
# sum_children n_r r^tau = sum_arrays k_w^tau), with W and targets passed in. Pure-target {73: W} is the hub's root.
from decimal import Decimal as D, getcontext
getcontext().prec = 60
_ln = {}
def _l(x):
    if x not in _ln: _ln[x] = D(x).ln()
    return _ln[x]
def root(children, targets, lo=D('0.97'), hi=D(1), iters=200):
    ch = [(D(n), _l(r)) for r, n in children.items() if n]
    tg = [(D(n), _l(k)) for k, n in targets.items() if n]
    for _ in range(iters):
        x = (lo + hi) / 2
        lhs = sum(n * (l * x).exp() for n, l in ch)
        rhs = sum(n * (l * x).exp() for n, l in tg)
        if lhs > rhs: lo = x
        else: hi = x
    return 1 - (lo + hi) / 2
def add(d, r, n):
    d = dict(d)
    if n: d[r] = d.get(r, 0) + n
    return d
def sep(children, targets):
    return max(r for r, n in children.items() if n) < min(k for k, n in targets.items() if n)
