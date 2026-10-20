#!/usr/bin/env python3
"""vrc.py (e3-rmbuild): exact rank-cover test (ex2-families Theorem 17's inequalities) for VIRTUAL slice systems
(no-phase game with a cyclic address group, Lemma C/V).  Fixed virtual multiset c_w, virtual arrays W n with parity
n_odd(Q).  For every form G: rank E[D_G, D_G] = slack(G) := |D_G| - W n_odd(G)  (Theorem 10 for projector targets).
For every form F and every partition of supp(F) x supp(F) into rectangles R x C (here: all partitions of supp F into
at most 3 blocks, pieces = block x block):  slack(F) <= sum over pieces of ub(R x C),
ub(R x C) = min(c(R), c(C), min{slack(G) : supp G contains R u C}).  A violation proves the target impossible.
Usage: vrc.py m n power W dirs."""
import sys, itertools
exec(open("vclique.py").read().split("m, n, power, W =")[0].split('"""', 2)[2])   # reuse the tf helpers
m, n, power, W = (int(x) for x in sys.argv[1:5])
dirs = [int(d) for d in sys.argv[5].split(",")]
N = 1 << m
sig = tf.singer_table(m, power) if n > 1 else list(range(N))
inv = [0] * N
for v in range(N):
    inv[sig[v]] = v
sk = [list(range(N))]
for k in range(1, n):
    sk.append([inv[v] for v in sk[-1]])
c = [0] * N
for z in dirs:
    for k in range(n):
        c[sk[k][z]] += 1
forms = [tuple(x // 2 for x in t) for t in tf.binary_forms(m)]
supp = []; slack = []
for Q in forms:
    S = frozenset(w for w in range(1, N) if Q[w])
    nodd = sum(1 for k in range(n) if sum(Q[sk[k][1 << j]] for j in range(m)) % 2)
    supp.append(S); slack.append(sum(c[w] for w in S) - W * nodd)
assert min(slack) >= 0, "deletion fails"
best_ub = {}
def ub(U):
    U = frozenset(U)
    if U not in best_ub:
        best_ub[U] = min((slack[g] for g in range(len(forms)) if U <= supp[g]), default=10 ** 9)
    return best_ub[U]
def piece(R, C):
    return min(sum(c[w] for w in R), sum(c[w] for w in C), ub(set(R) | set(C)))
def partitions(P, kmax):
    P = list(P)
    for k in range(1, kmax + 1):
        for lab in itertools.product(range(k), repeat=len(P) - 1):
            lab = (0,) + lab
            if max(lab) != k - 1 or any(lab[i] > max(lab[:i]) + 1 for i in range(1, len(lab))):
                continue
            yield [[P[i] for i in range(len(P)) if lab[i] == b] for b in range(k)]
viol = []; checked = 0
for f, S in enumerate(supp):
    if not S or slack[f] == 0:
        continue
    P = sorted(S)
    kmax = 3 if len(P) <= 8 else 2
    bestcover = None
    for part in itertools.chain(partitions(P, kmax), [[[x] for x in P]]):
        tot = sum(piece(R, C) for R in part for C in part)
        if bestcover is None or tot < bestcover:
            bestcover = tot
    checked += 1
    if bestcover < slack[f]:
        viol.append((f, len(P), slack[f], bestcover))
print("m=%d n=%d W=%d M=%d: forms checked %d, violations %d %s" % (m, n, W, len(dirs), checked, len(viol), viol[:4]))
print("VERDICT:", "EXCLUDED (rank cover)" if viol else "not excluded by rank cover with <= 3 blocks")
