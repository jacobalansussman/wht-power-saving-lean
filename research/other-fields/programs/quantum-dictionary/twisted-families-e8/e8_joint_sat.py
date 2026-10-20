#!/usr/bin/env python3
"""e8_joint_sat.py (t-family) -- NON-UNIFORM twists on the hub's complex P112 frame decoder: can any total cost less
than 8 when every total carries its own twists?

Relaxation (generous; UNSAT is an obstruction, SAT is only a floor-level possibility):
  total a (a = 0..7) reads and serves A_a = lines with coordinate a != 0, coefficients c_a[T] = L[T][a] (+-1).
  In total a, N_a[T,S] = omega(p'_aT, q'_aS) in F_2 (twisted pairing; t-sandwich), N_a[S,S] = 1, and the total
  costs >= rank_F2 N_a[A_a, A_a]: encoded as N_a = P_a^T Q_a with P_a, Q_a in F_2^{k_a x 56}.
  For a pair S != T the totals C = supp S cap supp T both carry x_T into Y_S.  They split into twist classes
  (equal relative twist); x_T terms of one class add as scalars, so a class whose sum of c_a[S] c_a[T] is non-zero
  must be cancelled directly by X_T (N_a[T,S] = 0 on that class, and N_a[S,T] = 0 by the same argument for x_S in
  Y_T); a zero-sum class needs nothing but has one twist, so N_a[T,S] and N_a[S,T] are constant on it.
  Allowed assignments of (N_a[T,S], N_a[S,T])_{a in C} = those for which SOME partition of C works.
  Ignored (more generous than reality): transitivity of twist classes across pairs, frame intervals, timing.
usage: python3 -I e8_joint_sat.py K0 [Kother] [totals, e.g. 0,1,2; default all 8] [--loose]
       python3 -I e8_joint_sat.py min BUDGET [--loose]   (every total rank 7 + e_a, sum e_a <= BUDGET)   (rank bound for total 0, rank bound for the others; default 7 8)"""
import os, sys, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '../lib')))
from copyfloor import e8_lines
from pysat.solvers import Cadical153
from pysat.formula import IDPool

L, G = e8_lines()
live = list(range(8, 120))
MODE_MIN = len(sys.argv) > 1 and sys.argv[1] == 'min'
if MODE_MIN:
    BUDGET = int(sys.argv[2]); sys.argv = [sys.argv[0], str(7 + BUDGET), str(7 + BUDGET)] + sys.argv[3:]
k0 = int(sys.argv[1]) if len(sys.argv) > 1 else 7
ko = int(sys.argv[2]) if len(sys.argv) > 2 else 8
LOOSE = '--loose' in sys.argv
if LOOSE: sys.argv.remove('--loose')
TOT = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else list(range(8))
A = {a: [S for S in live if L[S][a] != 0] for a in range(8)}
supp = {S: [a for a in range(8) if L[S][a] != 0] for S in live}

def partitions(s):
    s = list(s)
    if not s: yield []; return
    first, rest = s[0], s[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i+1:]
        yield [[first]] + p

def allowed(S, T, C):
    """set of allowed tuples (n_TS[a], n_ST[a]) for a in C (as a 2|C|-bit tuple)."""
    prod = {a: L[S][a] * L[T][a] for a in C}
    ok = set()
    for p in partitions(C):
        free = [cls for cls in p if sum(prod[a] for a in cls) == 0]
        # zero-sum classes: constant values (two bits each); others: zeros
        if LOOSE:     # --loose: values inside a zero-sum class need not agree (drops "N is a function of the twist")
            fr = [a for c_ in free for a in c_]
            for bits in itertools.product((0, 1), repeat=2 * len(fr)):
                val = {a: (0, 0) for a in C}
                for i, a in enumerate(fr): val[a] = (bits[2*i], bits[2*i+1])
                ok.add(tuple(x for a in C for x in val[a]))
            continue
        for bits in itertools.product((0, 1), repeat=2 * len(free)):
            val = {a: (0, 0) for a in C}
            for i, cls in enumerate(free):
                for a in cls: val[a] = (bits[2*i], bits[2*i+1])
            ok.add(tuple(x for a in C for x in val[a]))
    return ok

pool = IDPool(); cls = []
def AND(x, y):
    z = pool.id(('and', x, y)); cls.extend([[-z, x], [-z, y], [z, -x, -y]]); return z
def XOR_eq(lits, out):
    acc = lits[0]
    for l in lits[1:]:
        z = pool.id(('x', acc, l)); cls.extend([[-z, acc, l], [-z, -acc, -l], [z, -acc, l], [z, acc, -l]]); acc = z
    cls.extend([[-out, acc], [out, -acc]])

t0 = time.time()
Nv = {}
IND = []
def forced2(S, T):      # forced known in every class structure: non-orthogonal and sharing exactly two coordinates
    return G[S][T] != 0 and len([a for a in supp[S] if a in supp[T]]) == 2
def clique7(a):
    """the 7-clique K_sigma (sigma = all +): in each of the 7 supports through a, the line whose signs agree with
    sigma on every coordinate; pairwise non-orthogonal and sharing exactly two coordinates."""
    K = []
    for S in A[a]:
        sg = L[S][a]
        if all(L[S][c] * sg in (0, 1) for c in range(8)):
            K.append(S)
    assert len(K) == 7 and all(forced2(S, T) for S, T in itertools.combinations(K, 2)), K
    return K
# rank <= k parametrisation with symmetry broken by the clique K_a (N_a[K,K] = I_7, see the docstring of
# e8_minrank7.py; for k = 8 one extra rank-one term gamma_U delta_V with gamma = delta = 0 on K):
#   N_a[U,V] = sum_{S in K} N_a[U,S] N_a[S,V]  (+ gamma_U delta_V)
for a in TOT:
    k = k0 if a == 0 else ko
    K = clique7(a)
    for T in A[a]:
        for S in A[a]:
            Nv[(a, T, S)] = pool.id(('N', a, T, S))
    for S in K:
        for T in K:
            cls.append([Nv[(a, T, S)]] if S == T else [-Nv[(a, T, S)]])
    E = k - 7
    gam = [{U: pool.id(('g', a, t, U)) for U in A[a]} for t in range(E)]
    dlt = [{U: pool.id(('d', a, t, U)) for U in A[a]} for t in range(E)]
    for t in range(E):
        z = pool.id(('z', a, t)); IND.append(z)
        for U in K:
            cls.append([-gam[t][U]]); cls.append([-dlt[t][U]])
        for U in A[a]:
            cls.append([z, -gam[t][U]])           # unused term: gamma = 0
        if t:
            cls.append([-z, pool.id(('z', a, t - 1))])   # order the terms
    for U in A[a]:
        for V in A[a]:
            if U in K or V in K: continue
            lits = [AND(Nv[(a, U, S)], Nv[(a, S, V)]) for S in K] + [AND(gam[t][U], dlt[t][V]) for t in range(E)]
            XOR_eq(lits, Nv[(a, U, V)])
    for T in A[a]:
        cls.append([Nv[(a, T, T)]])
# symmetry breaking for total 0: a 7-clique of pairs that are forced known in every class structure
#   (non-orthogonal, sharing exactly two coordinates) gives N_0[K,K] = I; with k0 = 7 fix P_0[:,K] = Q_0[:,K] = I.
pairs = 0; forced = 0
for S, T in itertools.combinations(live, 2):
    Cfull = [a for a in supp[S] if a in supp[T]]
    C = [a for a in Cfull if a in TOT]
    if not C: continue
    okfull = allowed(S, T, Cfull); pairs += 1
    pos = [Cfull.index(a) for a in C]          # totals left out of TOT stay free: project the allowed set
    ok = {tuple(x for i in pos for x in (t[2*i], t[2*i+1])) for t in okfull}
    vars_ = [v for a in C for v in (Nv[(a, T, S)], Nv[(a, S, T)])]
    for bits in itertools.product((0, 1), repeat=len(vars_)):
        if bits not in ok:
            cls.append([(-v if b else v) for v, b in zip(vars_, bits)])
            forced += 1
print('k0 %d, others %d: %d vars, %d clauses (%d pair-forbids over %d pairs), built %.1fs'
      % (k0, ko, pool.top, len(cls), forced, pairs, time.time() - t0), flush=True)
if MODE_MIN:
    from pysat.card import CardEnc, EncType
    enc = CardEnc.atmost(lits=IND, bound=BUDGET, top_id=pool.top, encoding=EncType.seqcounter)
    cls.extend(enc.clauses)
    print('mode min: sum over totals of (rank - 7) <= %d, i.e. cst <= %d' % (BUDGET, 7 * len(TOT) + BUDGET), flush=True)
s = Cadical153(bootstrap_with=cls)
t1 = time.time(); r = s.solve()
print('SAT' if r else 'UNSAT', '(%.1fs)' % (time.time() - t1), flush=True)
if r:
    model = set(l for l in s.get_model() if l > 0)
    # verify the rank of each N_a over F_2 and the class conditions
    def rk(rows):
        basis = []
        for x in rows:
            for b in basis: x = min(x, x ^ b)
            if x: basis.append(x)
        return len(basis)
    for a in TOT:
        rows = [sum(1 << w for w, S in enumerate(A[a]) if Nv[(a, T, S)] in model) for T in A[a]]
        print('total', a, 'rank', rk(rows))
