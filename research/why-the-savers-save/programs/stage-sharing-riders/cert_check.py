"""r-riders / cert_check.py -- independent check of the 18-move bound (versions A and B), separate arithmetic.
Picks three DISJOINT flip pairs {z1,z1'}, {z2,z2'}, {z3,z3'} whose witnesses in lb_<v>.json have rank 6 (6 = number of
inputs, the maximum), then recomputes each rank from the END STATE OF THE REPLAYED PLAIN PROGRAM (replay.py, not the
target table of lb.py), evaluated with integers modulo p = 998244353 (p = 1 mod 4, i = a square root of -1 mod p).
A rank modulo p never exceeds the rank over Q(i); 6 columns cap it at 6; so mod-p rank 6 = rank 6 exactly.
Then n_z1 + n_z1' >= 6 three times over disjoint pairs gives total moves >= 18 for every word (Lemma A argument).
Also re-derives, for each witness pair, that every shift has equal value and exactly the pair flips.
usage: python3 -B cert_check.py"""
import sys, json, itertools
sys.dont_write_bytecode = True
import toy, replay
p = 998244353
I = next(pow(g, (p - 1) // 4, p) for g in range(2, 100) if pow(pow(g, (p - 1) // 4, p), 2, p) == p - 1)
assert I * I % p == p - 1
def chi(z): return [bin(y & z).count('1') & 1 for y in range(8)]
def value(q, m):   # i^(sum_y m_y q(y)), m over y = 1..7
    e = sum(m[y - 1] * q[y] for y in range(1, 8)) % 4
    return pow(I, e, p)
def rank_p(A):
    A = [r[:] for r in A]; rk = 0
    for c in range(len(A[0])):
        piv = next((i for i in range(rk, len(A)) if A[i][c] % p), None)
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]; inv = pow(A[rk][c], p - 2, p); A[rk] = [x * inv % p for x in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][c] % p:
                f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk
for v, prog in (('A', replay.prog_A_plain), ('B', replay.prog_B_plain)):
    M = replay.Machine('EQUAL', replay.LIVE); prog(M)
    assert not M.check(toy.targets(v)), 'replay end state differs from targets'
    rows = sorted(toy.targets(v))
    def mat(m):
        return [[sum(int(c) * value(q, m) for (inp, q), c in M.get(r).items() if inp == l) % p for l in toy.INPUTS] for r in rows]
    W = json.load(open('lb_%s.json' % v))['witnesses']
    six = {tuple(D): (m1, m2) for D, r, m1, m2 in W if len(D) == 2 and r == 6}
    trip = next(t for t in itertools.combinations(sorted(six), 3) if len(set(sum(t, ()))) == 6)
    print('[%s] disjoint pairs used: %s' % (v, list(trip)))
    for D in trip:
        m1, m2 = six[D]
        N1 = [sum(m1[y - 1] * chi(z)[y] for y in range(1, 8)) for z in range(1, 8)]
        N2 = [sum(m2[y - 1] * chi(z)[y] for y in range(1, 8)) for z in range(1, 8)]
        assert all((a - b) % 2 == 0 for a, b in zip(N1, N2)), 'a shift differs'
        flips = [z for z, a, b in zip(range(1, 8), N1, N2) if (a - b) % 4 == 2]
        assert flips == list(D), (flips, D)
        A, B = mat(m1), mat(m2)
        r = rank_p([[(a - b) % p for a, b in zip(ra, rb)] for ra, rb in zip(A, B)])
        print('     pair %s  points m=%s m\'=%s : shifts equal, flips exactly %s, rank of the end-state difference = %d'
              % (list(D), m1, m2, flips, r))
        assert r == 6
    print('[%s] => every word reaching this end state has >= 18 paid moves (6 per disjoint pair)' % v)
