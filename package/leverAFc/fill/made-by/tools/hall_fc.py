# fc-reopt (round wht45/fullcut): #3815's per-tick Hall LP (hub code, construction/hall_lp.py, loaded from a copy and
# re-parameterized in memory; the hub file is not changed) with the objective of a SEPARATED CUT of paid-rank cap C
# (C = 24: the full cut) instead of the hub's uniform objective. Generalises o-param's hall_any.py / o-model's
# mkpatch3815.py (half-cut head bonus) to:
#   * non-head bin with final remainder lf: cost (phi_C(lf) - lf)/a, phi_C(l) = floor(l/C) C^tau + (l mod C)^tau
#     (every paid block above C split, most unequal split = cheapest; tau = 1 - a);
#     --widen4: the final remainder is left - (number of width-4 roles) (#3945's offset step widens all 26 width-4
#     helpers to 5); patterns with left < n4 would overflow and get a prohibitive cost;
#   * EMPTY bin of a head class: cost (h_r - r)/a with h_r = (r - d)^tau + 73^tau - (73 - d)^tau, d = min(r, 72 - C)
#     (the head is cut by d axes, its array's target 73 -> 73 - d >= C + 1); head classes: pair class 4
#     (children 8/10, rank 48: the literal heads of the Lean build); with --ext also pair class 2 (children 4/6,
#     rank 32) and pair class 0 (children 0/2, rank 16), the other first-touch classes (o-model MODEL.md 4.4) [what-if:
#     not in our Lean];
#   * the same costs in the pricing step (reduced costs of every pattern).
# Up to the constant sum(lf)/a (fixed by the inventory) the objective is G(tau)/a, the separated-cut moment function
# at tau = 1 - a, so minimising it at a = b_k and re-rooting is Dinkelbach's method on the cut price.
import json, sys, types
import numpy as np

HEADCLASS = {4: 48, 2: 32, 0: 16}       # pair class -> rank (k = 9 r + i)


def phi(l, C, tau):
    if l <= 0: return 0.0
    q, m = divmod(int(l), C)
    return q * C ** tau + (m ** tau if m else 0.0)


def head_cost(r, C, tau):
    d = min(r, 72 - C)
    return ((r - d) ** tau if r - d else 0.0) + 73 ** tau - (73 - d) ** tau


def load(path, T, sup, local, C, ext, widen4):
    src = open(path).read()
    heads = {4} | ({0, 2} if ext else set())
    pre = ("PENALTY = 1e4\n"
           "import fc_cost as _fc\n"
           "FC_C = %d\nFC_HEADS = %r\nFC_WIDEN4 = %r\n" % (C, sorted(heads), widen4))
    rep = [("SUP = np.array([26, 2, 18, 53, 96, 213])\n", "SUP = np.array(%r)\n" % (sup,)),
           ("T = 360\n", "T = %d\n" % T),
           ("LOCAL = {1: 2440, 2: 725, 3: 224, 4: 102, 5: 20, 8: 8}\n", "LOCAL = %r\n" % (local,)),
           ("PENALTY = 1e4\n", pre),
           ("            costs.append(float(g_cost(left[j], self.a)))\n",
            "            costs.append(float(_fc.col_cost(P[j], int(left[j]), int(self.cap[k]), k % 9, self.a, FC_C, FC_HEADS, FC_WIDEN4)))\n"),
           ("            g = g_cost(left, self.a)\n",
            "            g = _fc.vec_cost(P, left, c, self.a, FC_C, FC_WIDEN4)\n"),
           ("            rc = g[None, :] - mu[ks, None] + pi[ks] @ P.T.astype(float)\n",
            "            rc = g[None, :] - mu[ks, None] + pi[ks] @ P.T.astype(float)\n"
            "            for _i in FC_HEADS:\n"
            "                if _fc.HEADCLASS[_i] == c:\n"
            "                    rc[(ks % 9 == _i)[:, None] & (left == c)[None, :]] += _fc.head_delta(c, self.a, FC_C)\n")]
    rep += [("            self.col_set.add((k, j))\n            self.cols.append((k, j))\n",
             "            for A, i, row in getattr(self, 'suffix', []):\n"
             "                ov = min(int(self.hi[k]), N - 1) - max(int(self.lo[k]) + 1, A) + 1\n"
             "                s = int(_fc.SLOTS[int(self.cap[k])][j, i])\n"
             "                if ov > 0 and s:\n"
             "                    ind.append(row)\n"
             "                    val.append(float(min(s, ov)))\n"
             "            self.col_set.add((k, j))\n            self.cols.append((k, j))\n"),
            ("            best = min(best, rc.min())\n",
             "            for A, i, row in getattr(self, 'suffix', []):\n"
             "                lam = self.dual[row]\n"
             "                if lam:\n"
             "                    ov = np.maximum(np.minimum(self.hi[ks], N - 1) - np.maximum(self.lo[ks] + 1, A) + 1, 0)\n"
             "                    rc -= lam * np.minimum(_fc.SLOTS[c][:, i][None, :], ov[:, None])\n"
             "            best = min(best, rc.min())\n")]
    for old, new in rep:
        assert src.count(old) == 1, old
        src = src.replace(old, new)
    m = types.ModuleType('hall_lp_fc')
    m.__file__ = path
    sys.modules['hall_lp'] = m
    exec(compile(src, 'hall_lp_fc', 'exec'), m.__dict__)
    import fc_cost
    fc_cost.SLOTS = {c: np.array([[sum(int(P[j, w]) for w in range(len(m.W)) if S >> w & 1) for S in m.SUBSETS]
                                  for j in range(len(P))]) for c, (P, _) in m.PATS.items()}

    def add_suffix(self, starts):
        """fc-reopt: for each suffix interval I = [A, N-1] and each set S of widths, the bins open during I can absorb at
        most min(slots of S in the bin, |I n window|) roles of I (one role per bin per call); the calls of I bring
        sum_{w in S} s_w |I| roles. Necessary for any placement (generalises the per-tick rows to intervals)."""
        from collections import defaultdict
        bycls = defaultdict(list)
        for n, (k, j) in enumerate(self.cols):
            bycls[k].append((j, self.first_y + n))
        self.suffix = getattr(self, 'suffix', [])
        lower, st, ind, val = [], [], [], []
        for A in starts:
            ov = np.maximum(np.minimum(self.hi, m.N - 1) - np.maximum(self.lo + 1, A) + 1, 0)
            ks = np.where(ov > 0)[0]
            for i in range(len(m.SUBSETS)):
                st.append(len(ind)); lower.append(float(m.SUBSET_SUP[i] * (m.N - A)))
                for k in ks:
                    for j, col in bycls[k]:
                        s = int(fc_cost.SLOTS[int(self.cap[k])][j, i])
                        if s:
                            ind.append(col); val.append(float(min(s, int(ov[k]))))
                self.suffix.append((int(A), i, self.nrows + len(lower) - 1))
        n = len(lower)
        self.h.addRows(n, np.array(lower), np.full(n, m.highspy.kHighsInf), len(ind), np.array(st, dtype=np.int32),
                       np.array(ind, dtype=np.int32), np.array(val))
        self.nrows += n
    m.HallLP.add_suffix = add_suffix
    return m


def summarize(H, lp, T, C, ext, widen4, final_local, W=None):
    """our cut price of the LP/MIP point: children = final local profile x 9T + paid pieces, targets per head."""
    tau0 = None
    y = lp.x[lp.first_y:]
    heads = {4} | ({0, 2} if ext else set())
    ch = {int(r): float(n) * 9 * T for r, n in final_local.items()}
    tg = {73: 896.0 * T}
    nhead = {}; frac = set(); lfh = {}
    for (k, j), v in zip(lp.cols, y):
        if v <= 1e-9: continue
        if abs(v - round(v)) > 1e-6: frac.add(k)
        c = int(lp.cap[k]); P = H.PATS[c][0][j]; left = int(H.PATS[c][1][j])
        i = k % 9
        if i in heads and left == c:
            d = min(c, 72 - C)
            if c - d: ch[c - d] = ch.get(c - d, 0) + v
            tg[73] -= v; tg[73 - d] = tg.get(73 - d, 0) + v
            nhead[c] = nhead.get(c, 0) + v
            continue
        lf = left - (int(P[0]) if widen4 else 0)
        lfh[lf] = lfh.get(lf, 0) + v
        if lf <= 0: continue
        q, m = divmod(lf, C)
        if q: ch[C] = ch.get(C, 0) + q * v
        if m: ch[m] = ch.get(m, 0) + v
    b = root(ch, tg)
    return dict(cut_price=b, heads=nhead, fractional_classes=len(frac), final_remainders={str(k): v for k, v in sorted(lfh.items())},
                children={str(k): v for k, v in sorted(ch.items())}, targets={str(k): v for k, v in sorted(tg.items())},
                slack=float(lp.slack().sum()), hall_rows=len(lp.hall), patterns=len(lp.cols))


def root(ch, tg):
    def G(t): return sum(n * r ** t for r, n in ch.items() if n) - sum(n * k ** t for k, n in tg.items() if n)
    lo, hi = 0.99, 0.99999
    for _ in range(200):
        mid = (lo + hi) / 2
        if G(mid) > 0: lo = mid
        else: hi = mid
    return 1 - (lo + hi) / 2
