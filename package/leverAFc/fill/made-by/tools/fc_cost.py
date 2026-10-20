# fc-reopt: column costs of the separated-cut objective for the patched hall_lp (see hall_fc.py). Units: hall_lp's
# g_a units, cost(l) = (moment(l) - l)/a, so the hub's PENALTY on unmet demand keeps its meaning.
import numpy as np

HEADCLASS = {4: 48, 2: 32, 0: 16}
BIG = 1e6


def phi(l, C, tau):
    if l <= 0: return 0.0
    q, m = divmod(int(l), C)
    return q * C ** tau + (m ** tau if m else 0.0)


def head_cost(r, C, tau):
    d = min(r, 72 - C)
    return ((r - d) ** tau if r - d else 0.0) + 73 ** tau - (73 - d) ** tau


import os as _os
LAM = float(_os.environ.get('FC_LAMBDA', '1'))     # 1: the cut objective; 0: the half-cut objective of design B
HC = 48                                           # the half cut leaves every paid rank at most 48: no split


def nonhead(lf, a, C):
    if lf < 0: return BIG
    full = (phi(lf, C, 1 - a) - lf) / a
    if LAM == 1: return full
    half = (phi(lf, HC, 1 - a) - lf) / a
    return LAM * full + (1 - LAM) * half


def col_cost(Prow, left, cap, cls, a, C, heads, widen4):
    if cls in heads and HEADCLASS.get(cls) == cap and left == cap:
        full = (head_cost(cap, C, 1 - a) - cap) / a
        if LAM == 1: return full
        return LAM * full + (1 - LAM) * (head_cost(cap, HC, 1 - a) - cap) / a
    if widen4 and int(Prow[0]) > 1: return BIG     # #3945's verify.py: each bin widens by at most one cell
    lf = left - (int(Prow[0]) if widen4 else 0)
    return nonhead(lf, a, C)


def vec_cost(P, left, cap, a, C, widen4):
    n4 = np.asarray(P)[:, 0].astype(np.int64) if widen4 else np.zeros(len(left), dtype=np.int64)
    lf = np.asarray(left, dtype=np.int64) - n4
    return np.array([BIG if (widen4 and k > 1) else nonhead(int(x), a, C) for x, k in zip(lf, n4)])


def head_delta(cap, a, C):
    tau = 1 - a
    full = (head_cost(cap, C, tau) - phi(cap, C, tau)) / a
    if LAM == 1: return full
    return LAM * full + (1 - LAM) * (head_cost(cap, HC, tau) - phi(cap, HC, tau)) / a
