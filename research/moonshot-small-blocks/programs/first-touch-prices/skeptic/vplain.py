# vplain.py -- consistency test of Lemma Q + the 105 table on GENUINE programs: the first touch of the plain programs
# is D1 x (letter = identity), so its 105-defect must fit the Theorem-M slack: plain R in 2 moves has slack 40 - 38 = 2,
# plain K in 3 has 60 - 49 = 11, plain 1+5 in 13 has 260 - 239 = 21.  Also: the identity letter's defect image, and
# the full over-demand OD of plain R and plain K in the 105 (exact, mod p), checked against the slack identity.
import sys, numpy as np
src = open(__file__.rsplit('/', 1)[0] + '/v105.py').read()
pre = src.split("    QhA = [mm(Q, mm(hW, W1(a), p), p) for a in Aall7]")[0]
pre = pre.replace("for p in PR:", "for p in PR[:1]:")
exec(pre)
I = np.eye(D, dtype=np.int64)
for nm in ('kernel', 'two-row'):
    print(f'identity letter: defect_{nm}(h) =', rank(mm(Q, mm(hW, E[nm], p), p), p))
# over-demand of the plain programs in the 105 at level 0: Kb = {v in N : pi C_t v = 0 for all t}
def OD(moves, T):
    Ccat = np.vstack([mm(Q, mm(M, NB, p), p) for M in moves])
    dKb = NB.shape[1] - rank(Ccat, p)
    cf = nullspace(mm(Q, mm(T, NB, p), p), p); dE = cf.shape[1]
    gains = []; cur = None
    for M in moves:
        blk = mm(Q, mm(M, NB, p), p); new = blk if cur is None else np.vstack([cur, blk])
        gains.append(rank(new, p) - (0 if cur is None else rank(cur, p))); cur = new
    return dE - dKb, gains
h2_ = W1(word_R7([ee(2)]))
od, g = OD([hW, mm(h2_, hW, p)], h2W); print('plain R in 2 (D1 then D2): OD', od, 'gains', g, 'slack identity 2*20 - 38 =', 40 - 38, '= sum(cap-gain)+OD =', sum(20 - x for x in g) + od)
h4_ = W1(word_R7([ee(3)]))
od, g = OD([hW, mm(h2_, hW, p), mm(h4_, mm(h2_, hW, p), p)], kW); print('plain K in 3: OD', od, 'gains', g, 'slack 60 - 49 = 11 = ', sum(20 - x for x in g) + od)
