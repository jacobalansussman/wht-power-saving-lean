"""r-riders / lb_prefix.py -- exact lower bound on the EXTRA paid moves once the 12 stage-i data climbs are fixed
(with or without v-s-wide's riders), wider game (hence also the equal-frame game).

Span form of Lemma A: for a formal vector v let d(v) = v(m) - v(m') in C^6 (its input coefficients evaluated at two
Walsh points m, m' with equal shift values).  Inputs have d = 0; mixing is linear; a shift multiplies d by the same
sign at m and m'; a move dir z with z NOT in the flip set D gives d(dir z c) = i^N_z d(c), already in the span; a move
with z in D adds at most one dimension.  So with S0 = C-span of d(outputs of the fixed climbs):
    #extra moves with direction in D  >=  dim(S0 + span d(target rows)) - dim(S0).
Interleaving does not matter: the climbs act on fixed contents (own data + the fixed riders).
Integer programme over the extra moves per direction as in lb.py; witnesses re-checked exactly (Fractions over Q(i)).
usage: python3 -B lb_prefix.py A|B plain|riders [pairs]"""
import sys, itertools, json
pass  # numpy/scipy: install with pip (the original run used a local copy)
sys.dont_write_bytecode = True
import numpy as np
import toy, lb, replay

def prefix_outputs(riders):
    M = replay.Machine('WIDE', replay.LIVE)
    outs = []
    orig_move = M.move
    def rec_move(r, z):
        orig_move(r, z); outs.append(dict(M.get(r)))
    M.move = rec_move
    replay.climbs(M, riders=riders)
    return outs

def ent_of(vec):
    assert all(c == int(c) for c in vec.values())
    return [[(int(c), q) for (inp, q), c in vec.items() if inp == l] for l in toy.INPUTS]

def run(v, mode):
    T = toy.targets(v); rows = sorted(T)
    tent = [ent_of(T[r]) for r in rows]
    pent = [ent_of(o) for o in prefix_outputs(mode == 'riders')]
    Mp = lb.all_points()
    TV = lb.eval_T(tent, Mp); PV = lb.eval_T(pent, Mp)      # (P, rows, 6), (P, 12, 6)
    idx = {tuple(m): k for k, m in enumerate(Mp.tolist())}
    chiM = np.array([[toy.chi(z)[y] for y in range(1, 8)] for z in toy.DIRS], dtype=np.int64)
    best = {}
    for d in itertools.product(range(4), repeat=7):
        nd = chiM @ np.array(d)
        if any(x % 2 for x in nd): continue
        D = frozenset(z for z, x in zip(toy.DIRS, nd % 4) if x == 2)
        if not D or (PAIRS and len(D) != 2): continue
        M2 = (Mp + np.array(d)) % 4
        k2 = np.array([idx[tuple(x)] for x in M2.tolist()])
        dP = PV - PV[k2]; dT = TV - TV[k2]
        r0 = np.linalg.matrix_rank(dP, tol=1e-7)
        r1 = np.linalg.matrix_rank(np.concatenate([dP, dT], axis=1), tol=1e-7)
        gain = r1 - r0; j = int(np.argmax(gain)); g = int(gain[j])
        if g > best.get(D, (-1,))[0]: best[D] = (g, tuple(Mp[j].tolist()), tuple(M2[j].tolist()))
    # exact re-check
    for D, (g, m1, m2) in best.items():
        def dm(ent):
            A = lb.exact_T(ent, m1); B = lb.exact_T(ent, m2)
            return [[(a[0] - b[0], a[1] - b[1]) for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]
        P0 = dm(pent); TT = dm(tent)
        eg = lb.exact_rank(P0 + TT) - lb.exact_rank(P0)
        assert eg == g, ('exact mismatch', sorted(D), g, eg)
    from scipy.optimize import milp, LinearConstraint, Bounds
    Ds = sorted(best, key=lambda D: (len(D), sorted(D)))
    A = np.array([[1 if z in D else 0 for z in toy.DIRS] for D in Ds]); lbv = np.array([best[D][0] for D in Ds])
    res = milp(c=np.ones(7), constraints=LinearConstraint(A, lbv, np.inf), integrality=np.ones(7), bounds=Bounds(0, np.inf))
    print('[%s, climbs %s] extra paid moves >= %d  (12 climbs + extra >= %d)   n_z = %s' % (v, mode, round(res.fun),
          12 + round(res.fun), dict(zip(toy.DIRS, [round(x) for x in res.x]))))
    print('    max gain per |D|: %s' % {s: max(best[D][0] for D in best if len(D) == s) for s in sorted({len(D) for D in best})})
    pairs = sorted([(sorted(D), best[D][0]) for D in best if len(D) == 2], key=lambda t: -t[1])
    print('    pair gains: %s' % pairs)
    json.dump({'version': v, 'mode': mode, 'extra_bound': round(res.fun),
               'witnesses': [[sorted(D), g, m1, m2] for D, (g, m1, m2) in best.items()]},
              open('lbp_%s_%s.json' % (v, mode), 'w'))

PAIRS = False
if __name__ == '__main__':
    PAIRS = len(sys.argv) > 3 and sys.argv[3] == 'pairs'   # only two-direction flip sets: fewer constraints, still a bound
    run(sys.argv[1], sys.argv[2])
