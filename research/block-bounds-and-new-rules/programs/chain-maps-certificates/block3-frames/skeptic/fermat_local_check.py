# v-xc3-cert3-compute: independent check of xc3-cert3's "Fermat-anchored lifts are capped at 2.287" (section 3.3).
# (i) orbits of the monomial group S_6 x (scaling a coordinate by alpha) on the 891 frames;
# (ii) for sampled other frames U: with the 405 Fermat lifts fixed, the necessary condition for U to be clean with a set J of
#      its Fermat neighbours (q_jk = Y_Vj cap Y_Vk 1-dim, all q_jk (j,k in J) of rank <= 3, and for each j the q_jk of rank <= 2),
#      tested EXACTLY over Q(omega) for every 4-subset J.  Also recheck: pairs of Fermat neighbours of U are at distance 2.
import itertools, sys, random
NSAMPLE = int(sys.argv[1]) if len(sys.argv) > 1 else 5
sys.argv = sys.argv[:1]
exec(open('lift891_check.py').read().split("def run(lift, name):")[0])   # reuse my own enumeration / Q(omega) code
def qnull(rows_cols):  # nullspace of the matrix whose COLUMNS are given vectors -> list of coefficient vectors
    vecs = rows_cols; n = len(vecs); m = len(vecs[0])
    M = [[(Fr(vecs[j][i][0]), Fr(vecs[j][i][1])) for j in range(n)] for i in range(m)]
    piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, m) if M[i][c] != Z0), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        inv = qinv(M[r][c]); M[r] = [qm(inv, v) for v in M[r]]
        for i in range(m):
            if i != r and M[i][c] != Z0:
                f = M[i][c]; M[i] = [qs(M[i][k], qm(f, M[r][k])) for k in range(n)]
        piv.append(c); r += 1
    free = [c for c in range(n) if c not in piv]
    basis = []
    for fc in free:
        v = [Z0] * n; v[fc] = (Fr(1), Fr(0))
        for i, pc in enumerate(piv): v[pc] = qs(Z0, M[i][fc])
        basis.append(v)
    return basis
def lin(coefs, vecs):
    out = [Z0] * 6
    for c, v in zip(coefs, vecs):
        for i in range(6): out[i] = (out[i][0] + qm(c, (Fr(v[i][0]), Fr(v[i][1])))[0], out[i][1] + qm(c, (Fr(v[i][0]), Fr(v[i][1])))[1])
    return out
# (i) orbits of the monomial group
def act_perm(U, perm):
    return frozenset(sum(((v >> (2 * perm[i])) & 3) << (2 * i) for i in range(6)) for v in U)
def act_scale(U, i):
    return frozenset(v ^ (((v >> (2 * i)) & 3) << (2 * i)) ^ (fm((v >> (2 * i)) & 3, 2) << (2 * i)) for v in U)
gens = [lambda U: act_perm(U, [1, 0, 2, 3, 4, 5]), lambda U: act_perm(U, [1, 2, 3, 4, 5, 0]), lambda U: act_scale(U, 0)]
seen = set(); orbits = []
for U in F:
    if U in seen: continue
    orb = {U}; st = [U]
    while st:
        X = st.pop()
        for g in gens:
            Y2 = g(X); assert Y2 in frames
            if Y2 not in orb: orb.add(Y2); st.append(Y2)
    seen |= orb; orbits.append(orb)
say(f'monomial group orbits on frames: sizes {sorted(len(o) for o in orbits)}; Fermat frames per orbit {[sum(1 for U in o if U in fermat) for o in orbits]}')
# Fermat lifts preserved by the generators (as subspaces over Q(omega)): perm obvious; scaling coordinate 0 by omega maps
# x_a = -z x_b to another Fermat plane -- checked: the set of 405 planes is closed (as frames) above; as lifts it is by formula.
# (ii) local test
liftF = {k: fermat[F[k]][0] for k in range(len(F)) if isF[k]}
others = [k for k in range(len(F)) if not isF[k]]
rng = random.Random(5)
sample = rng.sample(others, NSAMPLE)
for u in sample:
    J0 = [b for b in nbr[u] if isF[b]]
    dd = {}
    q = {}
    for j, k in itertools.combinations(J0, 2):
        d = dist(F[j], F[k]); dd[d] = dd.get(d, 0) + 1
        ns = qnull(liftF[j] + [[qs(Z0, (Fr(e[0]), Fr(e[1]))) for e in r] for r in liftF[k]])
        assert len(ns) == 1, len(ns)        # Y_Vj cap Y_Vk is exactly 1-dim (perfect family at distance 2)
        q[(j, k)] = q[(k, j)] = lin(ns[0][:3], liftF[j])
    best = 0; npass = 0
    for J in itertools.combinations(J0, 4):
        allq = [q[(j, k)] for j, k in itertools.combinations(J, 2)]
        ok = qrank(allq) <= 3 and all(qrank([q[(j, k)] for k in J if k != j]) <= 2 for j in J)
        if ok: npass += 1
    # largest J passing (3-subsets always pass? report)
    n3 = sum(1 for J in itertools.combinations(J0, 3)
             if qrank([q[(j, k)] for j, k in itertools.combinations(J, 2)]) <= 3 and all(qrank([q[(j, k)] for k in J if k != j]) <= 2 for j in J))
    say(f'other frame {u}: {len(J0)} Fermat nbrs; neighbour pairs by distance {dd}; 3-subsets passing the test {n3}/455; 4-subsets passing {npass}/1365 (exact)')
with open('out_fermat_local_check.txt', 'w') as fh: fh.write('\n'.join(out) + '\n')
