"""k1_sumzero.py (skeptic-2, wht44/ideas): which 4-in-15 survivors of rm2-phase A3(d) can exist WITHOUT scratch?

L-det (e7-rm, proved): with no scratch, sum over gates of z z^T = sum of the targets (mod 2). With each of the 15
directions of F_2^4 used once, the left side is 0, so the four targets must sum to 0.

Input: rm2-phase tile/surv_fe_alt_415_all.txt (1,612 classes up to GL(4,2); target indices into tile_alt's list of
the 448 nonsingular symmetric 4x4 matrices, alternating ones at index >= 420).
Independent check of the indexing: I rebuild the 448 matrices myself and compare with tile_alt's list.
"""
import sys, os, itertools
RM = os.environ.get('RM_PHASE_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                   'reed-muller-and-one-move', 'programs', 'phase-letters', 'ends-only-block4'))
RMD = os.environ.get('RM_PHASE_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                   'reed-muller-and-one-move', 'data', 'phase-letters', 'ends-only-block4'))
sys.path.insert(0, RM)
from tile_alt import target_matrices

m = 4
mats, nonalt = target_matrices(m)
print("targets:", len(mats), "non-alternating:", nonalt)

# independent rebuild: all symmetric 4x4 over F_2, nonsingular, split by alternating
def det2(M):
    M = [row[:] for row in M]; n = len(M); r = 0
    for c in range(n):
        p = next((i for i in range(r, n) if M[i][c]), None)
        if p is None: return 0
        M[r], M[p] = M[p], M[r]
        for i in range(n):
            if i != r and M[i][c]:
                M[i] = [a ^ b for a, b in zip(M[i], M[r])]
        r += 1
    return 1
mine = []
pairs = [(i, j) for i in range(m) for j in range(i, m)]
for bits in itertools.product((0, 1), repeat=len(pairs)):
    S = [[0] * m for _ in range(m)]
    for (i, j), b in zip(pairs, bits):
        S[i][j] = S[j][i] = b
    if det2(S): mine.append(S)
alt_mine = [S for S in mine if all(S[i][i] == 0 for i in range(m))]
key = lambda S: tuple(tuple(r) for r in S)
assert len(mine) == 448 and len(alt_mine) == 28
assert set(map(key, mine)) == set(map(key, mats))
assert all(all(mats[i][d][d] == 0 for d in range(m)) == (i >= nonalt) for i in range(len(mats)))
print("indexing check: 448 nonsingular, 28 alternating, alternating exactly at index >= %d: OK" % nonalt)

def add(A, B): return [[(a ^ b) for a, b in zip(r, s)] for r, s in zip(A, B)]
Z = [[0] * m for _ in range(m)]

rows = open(os.path.join(RMD, "surv_fe_alt_415_all.txt")).read().split()
by_nalt = {}; sumzero_by_nalt = {}; four_alt = []
for line in rows:
    dirs, tt = line.split("|")
    assert dirs == ",".join(str(i) for i in range(1, 16))
    t = [int(x) for x in tt.split(",")]
    na = sum(1 for x in t if x >= nonalt)
    S = Z
    for x in t: S = add(S, mats[x])
    sz = all(v == 0 for r in S for v in r)
    by_nalt[na] = by_nalt.get(na, 0) + 1
    if sz: sumzero_by_nalt[na] = sumzero_by_nalt.get(na, 0) + 1
    if na == 4: four_alt.append((t, sz, S))
print("classes by number of alternating targets:", dict(sorted(by_nalt.items())))
print("of which the four targets sum to 0 (L-det allows no-scratch):", dict(sorted(sumzero_by_nalt.items())))
print("total sum-zero classes:", sum(sumzero_by_nalt.values()), "of", len(rows))
print("the six four-alternating classes:")
for t, sz, S in four_alt:
    print("  targets", t, "sum zero" if sz else "sum = " + str(S))
