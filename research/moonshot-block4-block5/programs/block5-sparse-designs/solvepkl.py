"""solvepkl.py (w8-b5sparse): exact POINT tests on a pickled truncated d = 3 system (from alpha3's truncated mode):
for every V4 root choice, the layer-3 truncated system at random free parameters (two points), seed normalisation,
by sparse elimination with a dense FLINT finish.  Inconsistent at a random point = no solution for generic free
parameters on that component (special parameters not covered).  usage: solvepkl.py sys_d3_<seed>.pkl"""
import sys, pickle, random, itertools, time
import alpha3, v4torus
import modsolve as MS
import sparsemod

D = pickle.load(open(sys.argv[1], 'rb'))
nal = D['nal']; P2 = alpha3.P2
choices = []
if 'V' in D:
    V = D['V']
    for (a, b) in V['sols']:
        for tr in range(2):
            choices.append((('2col', a, b, tr), alpha3.alpha_point(V, nal, a, b, random.Random(100 * a + 10 * b + tr))))
else:
    T = D['T']
    ch = v4torus.root_choices(T)
    for ci, roots in enumerate(itertools.product(*ch)):
        for tr in range(2):
            r_ = random.Random(1000 + 10 * ci + tr)
            al = v4torus.alpha_from(T, roots, [r_.randrange(1, P2) for _ in range(T['N'] - T['k'])])
            al = [x if x is not None else r_.randrange(1, P2) for x in al]
            choices.append((('torus', ci, tr), al))
for (tag, al) in choices:
    out = []
    for l in (0, 1):
        cl, seedcols, cap = D['sys'][l]
        t = time.time()
        nr, nc, A = MS.build(cl, al, random.Random(7 + l), 'seed', seedcols, P2)
        ok, rk, st = sparsemod.consistent_sparse(A, {0: 1}, P2)
        out.append('out %d: %d x %d -> %s (rank %s, %s, %.0fs)' % (l + 1, nr, nc, {True: 'consistent', False: 'INCONSISTENT', None: 'gave up'}[ok], rk, st, time.time() - t))
    print(tag, ' | '.join(out), flush=True)
