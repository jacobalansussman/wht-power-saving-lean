"""v-r-swaps/indep.py -- skeptic's extra checks on r-swaps (13 Oct 2026).  usage: python3 -I indep.py
(1) A repair r-swaps did not try: carry the swapped content into stage 2 (swap back at local 0, free rank-0 route fix,
    as generous as B2f), then climb ONE rank unit to stage-2 local P_t, where X_t starts, and deliver Y_t += X_t with
    NO copy; then go home (via 0, or riding X_t to H and dropping the line).  Does any of this beat STD (168) or plain (160)?
(2) The parity lemma behind r-swaps' 'general reason', checked exhaustively on random subspaces over Q:
    for A <= B,  d(A,F) + d(F,B) - d(A,B) = 2[(rk F - dim F^B) + (rk A - dim A^F)], so a common frame F of the P-bank
    pair (P <= F <= I) and the 0-bank pair (0 <= F <= I-P) costs at least 2 extra unless P <= F <= I-P (impossible).
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.dont_write_bytecode = True
import swapdoor as sd
from fractions import Fraction as F

def carry_variant(home):
    M = sd.start_machine('carry-' + home)
    sd.forward_std(M, 0); M.reroute(0)
    sd.reverse_door(M, 1, 'B2')
    st1 = {k: dict(v) for k, v in M.cost.items()}
    M.reroute(1)
    for S_, T_ in sd.MATCH: M.swap('Y_' + S_, 'Y_' + T_, 'stage 2 swap back at local 0 (%s,%s)' % (S_, T_), back=True)
    # free rank-0 route fix (as generous as B2f): keep each term's own global frame, put it on rho(2,t)
    fixdist = []
    for t in sd.PORTS:
        R = M.reg['Y_' + t]; want = sd.key(sd.RHO[2, t])
        rs = {r for (u, r) in R['terms']}
        assert len(rs) == 1, rs
        r = [list(x) for x in next(iter(rs))]
        G = sd.mm(sd.mm(r, R['Lam']), sd.inv(r))
        Lfix = sd.mm(sd.mm(sd.RHOI[2, t], G), sd.RHO[2, t])
        R['terms'] = {(u, want): v for (u, rr), v in R['terms'].items()}; R['Lam'] = Lfix
        fixdist.append((sd.frame_dist(Lfix, sd.lift(2, sd.ZU)), sd.frame_dist(Lfix, sd.lift(2, sd.Pq(t))),
                        sd.nested(Lfix, sd.lift(2, sd.Pq(t)))))
    for t in sd.PORTS:
        M.move('Y_' + t, 'P:' + t, 'carry')                        # 1 rank unit (nested climb), checked by move()
        M.gate('Y_' + t, 'X_' + t, 1, 'stage 2 copy-free delivery %s at local P_%s' % (t, t))
        if home == 'via0':
            M.move('Y_' + t, '0', 'carry'); M.move('Y_' + t, 'K:' + t, 'climb')
        else:
            M.move('Y_' + t, 'H', 'climb'); M.move('Y_' + t, 'K:' + t, 'carry')
        M.move('X_' + t, 'H', 'climb')
    res = sd.final_check(M); tot = {k: v for k, v in M.cost.items() if not k.startswith('_')}
    return res, st1, tot, fixdist

def main():
    print('(1) carry the swapped content into stage 2 (free route fix, copy-free stage-2 delivery at local P_t)')
    for home in ('via0', 'rideH'):
        res, st1, tot, fd = carry_variant(home)
        ok = res['scalar_ok'] and res['routes_ok'] and res['frames_ok'] and res['kernel_ok']
        print('  %-6s passes all final checks: %s (off-route %d, junk %d); after the free route fix: distance to local 0 %s, '
              'to local P_t %s, nested in P_t %s' % (home, ok, res['off_route_terms'], res['junk_terms'],
              sorted({a for a, b, c in fd}), sorted({b for a, b, c in fd}), all(c for a, b, c in fd)))
        print('         paid moves whole word %s: rank mass %d (STD 168, plain 160, J 160)' % (tot, sd.mass(tot)))
    # (2) parity lemma on random rational subspaces of Q^6
    rnd = random.Random(7); n = 6; bad = 0; trials = 0; minexcess_meet = None
    def rand_sub(k, inside=None):
        if inside is None: vs = [[F(rnd.randint(-2, 2)) for _ in range(n)] for _ in range(k)]
        else:
            vs = []
            for _ in range(k):
                cs = [F(rnd.randint(-2, 2)) for _ in inside]
                vs.append([sum(c * w[i] for c, w in zip(cs, inside)) for i in range(n)])
        return vs
    def rk(vs): return sd.rank(vs) if vs else 0
    def meet(Avs, Bvs): return rk(Avs) + rk(Bvs) - rk(Avs + Bvs)
    def d(Avs, Bvs): return rk(Avs) + rk(Bvs) - 2 * meet(Avs, Bvs)
    for _ in range(400):
        B = rand_sub(rnd.randint(1, n)); A = rand_sub(rnd.randint(0, 3), B) if rk(B) else []
        Fv = rand_sub(rnd.randint(0, n))
        lhs = d(A, Fv) + d(Fv, B) - d(A, B)
        rhs = 2 * ((rk(Fv) - meet(Fv, B)) + (rk(A) - meet(A, Fv)))
        trials += 1; bad += lhs != rhs
    print('(2) parity lemma d(A,F)+d(F,B)-d(A,B) = 2[(rk F - dim F^B) + (rk A - dim A^F)] for A <= B: %d failures of %d' % (bad, trials))
    # the meeting cost on the tiny unit's U = Q^4 directly: min over all coordinate+label-spanned F of the pair excess
    q = [F(x) for x in sd.LABELS['w1']]; P = [q]; K = sd.perp_basis(q); Iv = [[F(int(i == j)) for j in range(4)] for i in range(4)]
    def d4(A, B):
        ra = sd.rank(A) if A else 0; rb = sd.rank(B) if B else 0; rab = sd.rank(A + B) if (A or B) else 0
        return ra + rb - 2 * (ra + rb - rab)
    cands = set(); vecs = [tuple(v) for v in Iv] + [tuple(q)] + [tuple(v) for v in K] + [(F(1), F(1), F(0), F(0)), (F(0), F(1), F(-1), F(1))]
    import itertools
    best = None
    for k in range(0, 5):
        for S in itertools.combinations(vecs, k):
            Fv = [list(v) for v in S]
            ex = (d4(P, Fv) + d4(Fv, Iv) - d4(P, Iv)) + (d4([], Fv) + d4(Fv, K) - d4([], K))
            best = ex if best is None else min(best, ex)
    print('    tiny unit, label w1: min pair excess of a common frame F over %d spanned candidates: %d (2 expected)' % (len(vecs), best))

if __name__ == '__main__':
    main()
