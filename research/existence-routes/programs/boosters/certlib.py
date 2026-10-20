# e3-boost certlib.py: the Theorem-R certificate profiles a relation of the job semiring must respect, and the screen.
#
# A relation rho (job model J, exist2/ex2-duality PROOF.md section 2) = tasks {k: w_k} served by calls {r: p_r}; net
# vector v_d = p_d - w_d.  Balance D(theta) = sum_d v_d d^theta.  A frame-symmetric certificate at block M (register of
# M rows) prices a rank-d job (paid block or target) at P_M(d) (Theorem R + Averaging Lemma; targets of any rank
# allowed), so every relation using at most M rows satisfies sum_d v_d P_M(d) >= 0.  Tensoring rho by [b] (rows ->
# b-row bundles, Lemma 2(b)) gives a relation on b*m rows with sizes b*d, so P_M(b d) applies whenever b*m <= M.
#
# Status of each profile: 'proved' (exact certificate with a written proof), 'whp' (exact integer data but f_j computed
# mod p, equal over Q with high probability), 'computed' (this folder, mod two primes), 'sampled' (averages from random
# pairs; NOT an exclusion), 'conj' (pattern verified for m <= 5 only).
from fractions import Fraction as F
from math import log

def A(M):  # label permutation module (w-bounds 5.4 A; v-ex2-duality: valid at every block)
    return [2 ** M - 2 ** (M - d) for d in range(1, M + 1)]

def B(M):  # A modulo equal frame sums (w-bounds 5.4 B); delta_B = delta_A - 1 computed for M = 2..5 (certB.py)
    return [2 ** M - 2 ** (M - d) - 1 for d in range(1, M + 1)]

PROFILES = []   # (name, M, [P(1..M)], status)
def add(name, M, P, status):
    assert len(P) == M, name
    PROFILES.append((name, M, [F(x) if not isinstance(x, float) else x for x in P], status))

for M in range(1, 65):
    add(f"A{M}", M, A(M), "proved")
for M in range(2, 6):
    add(f"B{M}", M, B(M), "computed" if M >= 4 else "proved")
for M in range(6, 65):
    add(f"B{M}", M, B(M), "conj")
# C_k: permutation module on the k-dimensional isotropic subspaces of the label space; Y_0 = those inside L_0.
# f(rank-d block) = G(m,k) - G(m-d,k) (Gaussian binomial, q = 2) [derived here: same proof as for A = C_1; at m = 3,
# C_1 = A (4,6,7) and C_2 = C (6,7,7) of w-bounds 5.4]; C_m = the frame count.
def G(n, k):
    if k < 0 or k > n: return 0
    num = den = 1
    for i in range(k):
        num *= 2 ** (n - i) - 1; den *= 2 ** (i + 1) - 1
    return num // den
def Ck(M, k):
    return [G(M, k) - G(M - d, k) for d in range(1, M + 1)]
for M in range(2, 41):
    for k in range(2, M):
        add(f"C{M}k{k}", M, Ck(M, k), "derived")

# Distance-set certificates on the frame graph (dualpolar.py): Y_U = span{e_W : dist(U, W) in S} [derived here].
import dualpolar as _dp
import itertools as _it
for M in range(2, 17):
    P_, ks_ = _dp.intersection_numbers(M)
    if M <= 8:
        Ss = [S for r in range(1, M + 1) for S in _it.combinations(range(M + 1), r)]
    else:
        Ss = [tuple(range(i, j + 1)) for i in range(M + 1) for j in range(i, M + 1) if (i, j) != (0, M)]
    for S in Ss:
        prof = _dp.dist_cert(M, S, P_, ks_)
        if max(prof) == 0: continue
        add(f"Dist{M}S{''.join(map(str, S)) if M <= 9 else '%d-%d' % (S[0], S[-1])}", M, prof, "derived")

# block 3, w-bounds 5.4 (smallreps.py, exact) and the hook / lifts
add("C3", 3, [6, 7, 7], "proved"); add("D3", 3, [5, 6, 6], "proved"); add("E3", 3, [12, 12, 15], "proved")
add("hook3", 3, [80, 144, 224], "proved")                                        # b56-compute pricecap.py (3,4)
add("liftE162", 3, [F(41, 35), F(1307, 630), 3], "proved")                       # w-bounds Theorem W2
add("liftg3c", 3, [F(1119, 945), F(7854, 3780), 3], "proved")                    # spectrum out_liftprof_m3.txt
# block 4
add("hook4", 4, [27456, 48576, 68320, 86592], "whp")                             # pricecap.py (4,8)
add("liftg4c", 4, [F(50595, 34425), F(785122, 321300), F(3512146, 1101600), 4], "proved")  # spectrum, all pairs exact
add("proj3to4", 4, [18480, 31920, 44352, 48720], "argued")                       # mixcap.py projection profile m=4
# block 5
add("hook5k3", 5, [6960, 10888, 12896, 13888, 14384], "whp")
add("hook5k4", 5, [64960, 103488, 124000, 133920, 138880], "whp")
add("proj3to5", 5, [18018000, 30597840, 41110272, 45002160, 47809440], "argued")
add("liftg5c", 5, [1.8758808834710163, 2.972255, 3.736826, 4.33294, 5.0], "sampled")
# block 6
add("hook6k3", 6, [60512, 92688, 108864, 116928, 120960, 122976], "whp")
add("lift6", 6, [2.374667, 3.572667, 4.458, 5.042667, 5.479333, 6.0], "sampled")

RIGOROUS = ("proved", "whp", "computed", "derived")

def D(v, t):
    return sum(float(c) * d ** t for d, c in v.items())

def window(v, lo=0.9, hi=0.9983, n=2000):
    """min of D on [lo, hi] (grid), argmin, and whether D(1) >= 0 (saves nothing alone)."""
    best = min(((D(v, lo + (hi - lo) * i / n), lo + (hi - lo) * i / n) for i in range(n + 1)))
    return best[0], best[1]

def roots(v, n=20000):
    """sign changes of D on [0, 1] (grid), as approximate roots."""
    out, prev = [], D(v, 0.0)
    for i in range(1, n + 1):
        t = i / n; cur = D(v, t)
        if (prev < 0) != (cur < 0): out.append(round(t, 5))
        prev = cur
    return out

def certify(v, m=None, statuses=RIGOROUS + ("argued", "sampled", "conj"), bmax=64):
    """All certificate tests for net vector v (dict rank -> weight) on a register of m rows (default max rank).
    Returns list of (name, b, status, margin) with margin = sum v_d P_M(b d); negative margin = excluded."""
    m = m or max(v)
    res = [("frame-count", 0, "proved", sum(F(c) if not isinstance(c, float) else c for c in v.values()))]
    for name, M, P, st in PROFILES:
        if st not in statuses: continue
        for b in range(1, bmax + 1):
            if b * m > M: break
            if (name[0] in "AB" or (name[0] == "C" and "k" in name) or name.startswith("Dist")) and b * m != M: continue   # A, B: use the smallest block b*m (B is strongest there)
            res.append((name, b, st, sum((c if isinstance(c, float) else F(c)) * P[b * d - 1] for d, c in v.items())))
    return res

def verdict(v, m=None, statuses=RIGOROUS + ("argued", "sampled", "conj")):
    res = certify(v, m, statuses)
    bad = [r for r in res if r[3] < 0]
    rig = [r for r in bad if r[2] in RIGOROUS + ("argued",)]
    return ("EXCLUDED" if rig else ("excluded-weakly" if bad else "PASSES")), bad, res

def mk(tasks, calls):
    v = {}
    for k, w in tasks.items(): v[k] = v.get(k, 0) - w
    for r, p in calls.items(): v[r] = v.get(r, 0) + p
    return {d: c for d, c in v.items() if c != 0}

# ---- the walk condition (derived here; e3-boost NOTE section 4) ----
# In the formal-content picture every term of a register is c [w] y_l, w = f_j P_j ... f_1 P_1 f_0 with P_i the paid
# steps that term passed through (each at most once, in time order) and f_i free letters (isometries of the frame
# graph fixing L_0).  A correct word must contain a term whose matrix is the target H_S itself (formal identity), so
# dist(L_0, w L_0) = |S| = k.  Along the word, x_i = dist(L_0, u_i L_0) moves by at most r_i (triangle inequality),
# so k lies in [max(0, 2 r_max - sum R'), min(sum R', M)] for the multiset R' of ranks on that term's path.
# Necessary condition: every target rank k is reachable from SOME sub-multiset R' of the calls.
from itertools import combinations as _comb
def walk_ok(v, M=None):
    M = M or max(v)
    calls = []
    for d, c in v.items():
        if c > 0: calls += [d] * int(c)
    tasks = sorted(d for d, c in v.items() if c < 0)
    reach = set()
    for r in range(1, len(calls) + 1):
        for R in set(_comb(sorted(calls), r)):
            lo, hi = max(0, 2 * max(R) - sum(R)), min(sum(R), M)
            reach.update(range(lo, hi + 1))
    return all(k in reach for k in tasks)

# ---- projections (b56-methods project.py): a frame family Y' at block M-1 gives Y_L = (+)_p Y'(pi_p L) at block M,
# abar_r = sum_s N_r(s) abar'_s (abar'_0 = 0).  Closed forms, matching project.py's exact counts at M = 4..7
# (proj/N5.txt .. N7.txt, mixcap.py N4..N6): N_r(r-2) = 2^(M-r)(2^r - 1)(2^(r-1) - 1),
# N_r(r-1) = 2^(M-r)(2^r - 1)(2^(r-1) + 2), N_r(r) = 4^M - 1 - the other two.  Status 'argued': rests on the Averaging
# Lemma at block M (Sp(2M,2) distance-transitive on frames; checked by matrices at M = 3), as b56-methods.
def Nproj(M):
    N = {}
    for r in range(1, M + 1):
        a = 2 ** (M - r) * (2 ** r - 1) * (2 ** (r - 1) - 1); b = 2 ** (M - r) * (2 ** r - 1) * (2 ** (r - 1) + 2)
        N[r] = {r - 2: a, r - 1: b, r: 4 ** M - 1 - a - b}
    return N
def project(prev):
    M = len(prev) + 1; N = Nproj(M); pv = [0] + list(prev)
    return [sum(c * pv[s] for s, c in N[r].items() if 0 <= s <= M - 1) for r in range(1, M + 1)]
assert project([80, 144, 224]) == [18480, 31920, 44352, 48720]          # = proj3to4 (mixcap.py)
_base = [(n, M, P, s) for (n, M, P, s) in list(PROFILES) if n in ("hook3", "liftE162", "liftg3c", "hook4", "liftg4c",
                                                                     "hook5k3", "hook5k4", "hook6k3", "E3", "B3", "C3", "D3")]
for (n, M, P, s) in _base:
    cur = list(P)
    for M2 in range(M + 1, 17):
        cur = project(cur)
        if M2 <= 8 or n in ("hook3", "liftE162", "liftg4c", "hook4"):
            add(f"P{n}to{M2}", M2, cur, "argued")
