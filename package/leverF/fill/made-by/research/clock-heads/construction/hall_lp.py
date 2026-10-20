"""Exact per-tick LP relaxation of the T360 owner fill for a given clock (diagnostic and
construction provenance; not run by verify.py). Comments are in French.

PL exact, instant par instant, du rangement T360 de #3734 pour une horloge donnée.

Classe = (réplique propriétaire, paire d'enfants) : 224 cases de même capacité et de même
fenêtre (lo, hi]. Chaque instant apporte le même lot de 408 objets ; pour une largeur w, le
transport instants → classes (intervalles) est réalisable si et seulement si, pour tout
intervalle d'instants [A, B], la demande des classes contenues dans [A, B] ne dépasse pas
l'offre s_w (B − A + 1) (Hall ; les totaux sont égaux). Ces lignes sont générées à la demande,
les motifs par classe aussi (tarification exacte sur tous les motifs de la capacité).

Relaxation : ni l'intégrité des motifs ni la règle « un objet par appel et par case » ne sont
imposées. L'optimum majore donc κ pour cette horloge ; il ne dit pas qu'un rangement entier
l'atteint.
"""
import gzip
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import highspy
import numpy as np

W = np.array([4, 5, 6, 7, 8, 9])
SUP = np.array([26, 2, 18, 53, 96, 213])
T = 360
N = 9 * T
LOCAL = {1: 2440, 2: 725, 3: 224, 4: 102, 5: 20, 8: 8}
BINS = 224
# (capacité, étape d'ouverture ou -1, étape de fermeture) des neuf paires, lues dans TEMPLATE
SPANS = [(16, -1, 1), (8, 1, 2), (32, -1, 3), (8, 3, 4), (48, -1, 5), (8, 5, 6), (48, 2, 8), (32, 4, 8), (16, 6, 8)]
PENALTY = 1e4
import math as _m
HEAD_SCALE = 1.0
def head_bonus(a):
    t = 1 - a
    return HEAD_SCALE * (48 ** t - 24 ** t - 73 ** t + 49 ** t) / a


def all_patterns(cap):
    out = []

    def rec(i, rest, cur):
        if i == len(W):
            out.append(tuple(cur))
            return
        for k in range(rest // W[i] + 1):
            rec(i + 1, rest - k * W[i], cur + [k])
    rec(0, cap, [])
    P = np.array(out, dtype=np.int64)
    return P, cap - P @ W


PATS = {c: all_patterns(c) for c in (8, 16, 32, 48)}
SUBSETS = [S for S in range(1, 64)]
SUBSET_SUP = np.array([sum(int(SUP[w]) for w in range(6) if S >> w & 1) for S in SUBSETS])
# INTER[c][j, i] = 1 si le motif j de capacité c a une place dans le sous-ensemble i
INTER = {c: np.array([[1.0 if any(P[j, w] and S >> w & 1 for w in range(6)) else 0.0 for S in SUBSETS]
                      for j in range(len(P))]) for c, (P, _) in PATS.items()}


def g_cost(left, a):
    left = np.asarray(left, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        v = (left ** (1 - a) - left) / a
    return np.where(left == 0, 0.0, v)


def F(profile, a):
    return sum(n * r ** (1 - a) for r, n in profile.items()) - 896 * T * 73 ** (1 - a)


def root(profile):
    lo, hi = 1e-4, 5e-3
    for _ in range(200):
        mid = (lo + hi) / 2
        if F(profile, mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo


def read_order(path):
    raw = Path(path).read_bytes()
    return json.loads(gzip.decompress(raw) if str(path).endswith('.gz') else raw)


def taus(order):
    tau = np.full((T, 9), -1, dtype=np.int64)
    for t, (s, r) in enumerate(order):
        tau[r, s] = t
    assert (tau >= 0).all() and (np.diff(tau, axis=1) > 0).all(), 'ordre des étapes par réplique'
    return tau


def classes(order):
    tau = taus(order)
    cap, lo, hi = [], [], []
    for r in range(T):
        for c, s0, s1 in SPANS:
            cap.append(c)
            lo.append(-1 if s0 < 0 else tau[r, s0])
            hi.append(tau[r, s1])
    return np.array(cap), np.array(lo), np.array(hi)


class HallLP:
    def __init__(self, order, a=0.0017304483602076, seeds=None, verbose=True):
        self.a = a
        self.cap, self.lo, self.hi = classes(order)
        self.K = len(self.cap)
        self.verbose = verbose
        h = highspy.Highs()
        h.setOptionValue('output_flag', False)
        h.setOptionValue('solver', 'simplex')
        self.h = h
        K = self.K
        # variables d[k, w] : 0 .. 6K-1 ; écarts u[w] : 6K .. 6K+5
        inf = highspy.kHighsInf
        h.addVars(6 * K, np.zeros(6 * K), np.full(6 * K, inf))
        h.addVars(6, np.zeros(6), np.full(6, inf))
        h.changeColsCost(6, np.arange(6 * K, 6 * K + 6, dtype=np.int32), np.full(6, PENALTY))
        # lignes : R1[k] somme des motifs = 224 ; R2[k,w] d - sum p y = 0 ; R3[w] total
        self.r1 = 0
        h.addRows(K, np.full(K, float(BINS)), np.full(K, float(BINS)), 0,
                  np.zeros(K, dtype=np.int32), np.array([], dtype=np.int32), np.array([]))
        self.r2 = K
        idx = np.arange(6 * K, dtype=np.int32)
        h.addRows(6 * K, np.zeros(6 * K), np.zeros(6 * K), 6 * K, idx, idx, np.ones(6 * K))
        self.r3 = 7 * K
        starts = np.arange(6, dtype=np.int32) * (K + 1)
        ind, val = [], []
        for w in range(6):
            ind += list(range(w, 6 * K, 6)) + [6 * K + w]
            val += [1.0] * (K + 1)
        tot = (SUP * N).astype(float)
        h.addRows(6, tot, tot, len(ind), starts, np.array(ind, dtype=np.int32), np.array(val))
        self.nrows = 7 * K + 6
        self.hall = []                    # (w, A, B)
        self.hall_rows = []               # indice de ligne de chaque condition intérieure
        self.nhall_flushed = 0
        self.hall_set = set()
        self.distinct = []                # (t, i sous-ensemble, ligne)
        self.cols = []                    # (k, index du motif)
        self.col_set = set()
        self.ncols = 6 * K + 6
        # motifs de départ : vide + motifs à reste nul ou idéal, pour chaque capacité
        if seeds is None:
            seeds = {}
            for c, (P, left) in PATS.items():
                keep = np.where((P.sum(1) == 0) | (left <= 1) | ((c == 32) & (left == 5) & (P[:, 5] == 3))
                                | ((c == 48) & (left == 3) & (P[:, 5] == 5)))[0]
                seeds[c] = keep
        self.add_cumulative()
        new = [(k, int(j)) for k in range(K) for j in seeds[self.cap[k]]]
        self.add_columns(new)

    def add_cumulative(self):
        """Lignes de Hall préfixe [0, B] et suffixe [A, N-1], creuses grâce à des cumuls.

        C_w(b_i) = C_w(b_{i-1}) + somme des d[k, w] des classes fermées en b_i, borné par s_w (b_i + 1) ;
        S_w(a_i) = S_w(a_{i+1}) + somme des d[k, w] des classes ouvertes en a_i, borné par s_w (N - a_i).
        """
        h = self.h
        lo1 = self.lo + 1
        ends, starts = np.unique(self.hi), np.unique(lo1)
        for w in range(6):
            for keys, kind in ((ends, 'C'), (starts, 'S')):
                n = len(keys)
                first = self.ncols
                ub = SUP[w] * (keys + 1) if kind == 'C' else SUP[w] * (N - keys)
                h.addVars(n, np.zeros(n), ub.astype(float))
                self.ncols += n
                rs, ind, val = [], [], []
                for i, b in enumerate(keys):
                    rs.append(len(ind))
                    ind.append(first + i)
                    val.append(1.0)
                    prev = i - 1 if kind == 'C' else i + 1
                    if 0 <= prev < n:
                        ind.append(first + prev)
                        val.append(-1.0)
                    ks = np.where((self.hi if kind == 'C' else lo1) == b)[0]
                    ind.extend((6 * ks + w).tolist())
                    val.extend([-1.0] * len(ks))
                h.addRows(n, np.zeros(n), np.zeros(n), len(ind), np.array(rs, dtype=np.int32),
                          np.array(ind, dtype=np.int32), np.array(val))
                self.nrows += n
                for b in keys:
                    self.hall_set.add((w, 0, int(b)) if kind == 'C' else (w, int(b), N - 1))
        self.first_y = self.ncols

    def add_columns(self, new):
        new = [x for x in new if x not in self.col_set]
        if not new:
            return 0
        costs, starts, ind, val = [], [], [], []
        for k, j in new:
            P, left = PATS[self.cap[k]]
            costs.append(float(g_cost(left[j], self.a)) - (head_bonus(self.a) if (k % 9 == 4 and int(left[j]) == 48) else 0.0))
            starts.append(len(ind))
            ind.append(self.r1 + k)
            val.append(1.0)
            for w in range(6):
                if P[j, w]:
                    ind.append(self.r2 + 6 * k + w)
                    val.append(-float(P[j, w]))
            for t, i, row in self.distinct:
                if self.lo[k] < t <= self.hi[k] and INTER[self.cap[k]][j, i]:
                    ind.append(row)
                    val.append(1.0)
            self.col_set.add((k, j))
            self.cols.append((k, j))
        n = len(new)
        self.h.addCols(n, np.array(costs), np.zeros(n), np.full(n, highspy.kHighsInf), len(ind),
                       np.array(starts, dtype=np.int32), np.array(ind, dtype=np.int32), np.array(val))
        self.ncols += n
        return n

    def add_distinct(self, ticks):
        """À chaque instant t : pour tout sous-ensemble S de largeurs, au moins sum_{w in S} s_w cases
        ouvertes à t ont une place dans S (couplage par appel ; nécessaire, pas suffisant)."""
        bycls = defaultdict(list)
        for n, (k, j) in enumerate(self.cols):
            bycls[k].append((j, self.first_y + n))
        lower, starts, ind, val = [], [], [], []
        for t in ticks:
            ks = np.where((self.lo < t) & (t <= self.hi))[0]
            for i in range(len(SUBSETS)):
                starts.append(len(ind))
                lower.append(float(SUBSET_SUP[i]))
                for k in ks:
                    for j, col in bycls[k]:
                        if INTER[self.cap[k]][j, i]:
                            ind.append(col)
                            val.append(1.0)
                self.distinct.append((t, i, self.nrows + len(lower) - 1))
        n = len(lower)
        self.h.addRows(n, np.array(lower), np.full(n, highspy.kHighsInf), len(ind),
                       np.array(starts, dtype=np.int32), np.array(ind, dtype=np.int32), np.array(val))
        self.nrows += n

    def queue_hall(self, w, A, B):
        if (w, A, B) in self.hall_set or B < A:
            return
        self.hall_set.add((w, A, B))
        self.pending = getattr(self, 'pending', [])
        self.pending.append((w, A, B))

    def flush_hall(self):
        rows = getattr(self, 'pending', [])
        self.pending = []
        if not rows:
            return 0
        lo1 = self.lo + 1
        lower, upper, starts, ind, val = [], [], [], [], []
        for w, A, B in rows:
            ks = np.where((lo1 >= A) & (self.hi <= B))[0]
            starts.append(len(ind))
            ind.extend((6 * ks + w).tolist())
            val.extend([1.0] * len(ks))
            lower.append(-highspy.kHighsInf)
            upper.append(float(SUP[w] * (B - A + 1)))
            self.hall.append((w, A, B))
            self.hall_rows.append(self.nrows + len(self.hall_rows) - self.nhall_flushed)
        n = len(rows)
        self.h.addRows(n, np.array(lower), np.array(upper), len(ind), np.array(starts, dtype=np.int32),
                       np.array(ind, dtype=np.int32), np.array(val))
        self.nrows += n
        self.nhall_flushed = len(self.hall_rows)
        return n

    def solve(self):
        self.h.run()
        st = self.h.getModelStatus()
        assert st == highspy.HighsModelStatus.kOptimal, self.h.modelStatusToString(st)
        sol = self.h.getSolution()
        self.x = np.array(sol.col_value)
        self.dual = np.array(sol.row_dual)
        self.obj = self.h.getInfo().objective_function_value
        return self.obj

    def demand(self):
        return self.x[:6 * self.K].reshape(self.K, 6)

    def separate(self, per_width=150, tol=1e-6):
        d = self.demand()
        lo1, hi = self.lo + 1, self.hi
        added = 0
        Bidx = np.arange(N)
        for w in range(6):
            M = np.zeros((N, N))
            np.add.at(M, (lo1, hi), d[:, w])
            D = np.cumsum(np.cumsum(M[::-1], axis=0)[::-1], axis=1)
            V = D - SUP[w] * (Bidx[None, :] - Bidx[:, None] + 1)
            V[np.tril_indices(N, -1)] = -np.inf
            bestB = V.argmax(1)
            bestV = V[np.arange(N), bestB]
            order = np.argsort(-bestV)
            cnt = 0
            for A in order:
                if bestV[A] <= tol or cnt >= per_width:
                    break
                if (w, int(A), int(bestB[A])) not in self.hall_set:
                    self.queue_hall(w, int(A), int(bestB[A]))
                    cnt += 1
            added += cnt
            del M, D, V
        self.flush_hall()
        return added

    def price(self, per_class=3, tol=1e-7):
        mu = self.dual[self.r1:self.r1 + self.K]
        pi = self.dual[self.r2:self.r2 + 6 * self.K].reshape(self.K, 6)
        new = []
        best = 0.0
        for c, (P, left) in PATS.items():
            ks = np.where(self.cap == c)[0]
            g = g_cost(left, self.a)
            rc = g[None, :] - mu[ks, None] + pi[ks] @ P.T.astype(float)
            if c == 48:
                rc[(ks % 9 == 4)[:, None] & (left == 48)[None, :]] -= head_bonus(self.a)
            for t, i, row in self.distinct:
                lam = self.dual[row]
                if lam:
                    open_ = (self.lo[ks] < t) & (t <= self.hi[ks])
                    rc[open_] -= lam * INTER[c][:, i][None, :]
            best = min(best, rc.min())
            for i, k in enumerate(ks):
                js = np.argsort(rc[i])[:per_class]
                for j in js:
                    if rc[i, j] < -tol and (k, int(j)) not in self.col_set:
                        new.append((int(k), int(j)))
        return self.add_columns(new), best

    def profile(self):
        prof = {r: float(N * n) for r, n in LOCAL.items()}
        for (k, j), v in zip(self.cols, self.x[self.first_y:]):
            if v > 1e-12:
                l = int(PATS[self.cap[k]][1][j])
                if l:
                    prof[l] = prof.get(l, 0.0) + v
        return prof

    def slack(self):
        return self.x[6 * self.K:6 * self.K + 6]

    def run(self, max_rounds=200):
        t0 = time.time()
        for it in range(max_rounds):
            self.solve()
            nr = self.separate()
            nc, best = self.price()
            prof = self.profile()
            if self.verbose:
                print(f'  tour {it:3d} obj {self.obj:.6f} écart {self.slack().sum():.3g} +{nr} lignes +{nc} motifs '
                      f'(rc min {best:.2e}) lignes Hall {len(self.hall)} motifs {len(self.cols)} '
                      f'racine {root(prof):.16f} [{time.time() - t0:.0f}s]', flush=True)
            if nr == 0 and nc == 0:
                break
        return root(self.profile())


IDEAL = 0.00173057843726628
CURRENT = 0.00173044836020758


def report(lp, label):
    prof = lp.profile()
    r = root(prof)
    idle = {k: round(v - N * LOCAL.get(k, 0), 3) for k, v in sorted(prof.items()) if abs(v - N * LOCAL.get(k, 0)) > 1e-6}
    gap_units = F(prof, IDEAL)
    print(f'{label}: racine PL {r:.16f}  idéal {IDEAL:.16f}  actuel {CURRENT:.16f}')
    print(f'   écart à l\'idéal {IDEAL - r:.3e} ; gain possible sur l\'actuel {r - CURRENT:+.3e}')
    print(f'   restes {idle}  (F à la racine idéale {gap_units:.4f})')
    return r


if __name__ == '__main__':
    order = read_order(sys.argv[1])
    lp = HallLP(order)
    lp.run()
    report(lp, sys.argv[1])
