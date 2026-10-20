# independent recheck of d-g-crossstage section (3)/(4): real subspace tests in the agent's own toy model
import sys, os, itertools
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-floor-checker')))
import copyfloor as cf
def rank(vs):
    basis = {}
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h in basis: v ^= basis[h]
            else: basis[h] = v; break
    return len(basis)
def dot(a, b): return bin(a & b).count('1') & 1
def perp9(vs):
    allv = [x for x in range(1, 512) if all(dot(x, c) == 0 for c in vs)]
    b = []
    for x in allv:
        if rank(b + [x]) > len(b): b.append(x)
    return b
L, G = cf.e8_lines()
B0 = [[(1 + G[i][j]) & 1 for j in range(120)] for i in range(120)]
U = cf.f2_labels(B0)
u = [sum(U[S][c] << c for c in range(len(U[S]))) for S in range(120)]
live = list(range(8, 120))
print('norms of live labels:', set(dot(u[S], u[S]) for S in live), ' distinct labels', len(set(u[S] for S in live)))
orth = [sum(1 for T in live if T != S and dot(u[S], u[T]) == 0) for S in live]
print('orthogonal neighbours per port:', set(orth), ' total off-diagonal orthogonal pairs', sum(orth))
# toy ambient F2^17 = U (bits 0..8) (+) B (bits 9..16), everything relative to E_i
Ub = [1 << k for k in range(9)]; Bb = [1 << (9 + k) for k in range(8)]
def contains(big, small): return rank(big + small) == rank(big)
def common_frame(lo1, hi1, lo2, hi2):
    # exists G with lo1<=G<=hi1, lo2<=G<=hi2  iff  lo1+lo2 <= hi1 cap hi2  iff lo1+lo2 <= hi1 and lo1+lo2 <= hi2
    lo = lo1 + lo2
    return contains(hi1, lo) and contains(hi2, lo)
cnt_fr = cnt_rf = 0
for S in live:
    for T in live:
        aTp = perp9([u[T]]); aSp = perp9([u[S]])
        # fwd->rev: target P-bank S in [U, U+B]; source 0-bank T in [aT^perp, aT^perp+B]
        if common_frame(Ub, Ub + Bb, aTp, aTp + Bb): cnt_fr += 1
        # rev->fwd: target 0-bank S in [aS^perp, aS^perp+B]; source P-bank T in [U, U+B]
        if common_frame(aSp, aSp + Bb, Ub, Ub + Bb): cnt_rf += 1
print('OLD identification, real subspace test: fwd->rev %d, rev->fwd %d of %d' % (cnt_fr, cnt_rf, 112 * 112))
# does the obstruction depend on N(a,a)=1?  try an isotropic label (even weight) in the same test
iso = 0b11
aTp = perp9([iso]); print('isotropic control a=%s: a in a^perp? %s; common frame? %s' % (bin(iso), contains(aTp, [iso]), common_frame(Ub, Ub + Bb, aTp, aTp + Bb)))
# NEW identification: target P-bank S in [<a_S>, U]; source 0-bank T in [0, a_T^perp]
cnt_new = sum(1 for S in live for T in live if common_frame([u[S]], Ub, [], perp9([u[T]])))
cnt_new_diag = sum(1 for S in live if common_frame([u[S]], Ub, [], perp9([u[S]])))
print('NEW identification, real subspace test: %d pairs (diagonal %d)' % (cnt_new, cnt_new_diag))
# section (4): totals A_a for a in 0..7; dim of read span, and dim of perp
for a in range(8):
    A = [S for S in live if G[S][a] != 0]
    P = [u[S] for S in A]
    print('total', a, 'served/read', len(A), 'w_in', rank(P), 'dim P^perp', len(perp9(P)))
