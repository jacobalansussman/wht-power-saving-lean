"""localize.py (w8-b5sparse): find a SMALL valid subsystem of a pickled truncated d = 3 system that is already
inconsistent at random free parameters: rows within BFS distance k of the seed column (with every column touching
them; any set of complete rows is a valid test), seed normalisation.  Prints the smallest k, the subsystem size, and
which free-parameter monomials occur in it.  usage: localize.py sys.pkl out(1|2) rootchoice kmax"""
import sys, pickle, random, itertools
import alpha3, v4torus
import modsolve as MS
import sparsemod

D = pickle.load(open(sys.argv[1], 'rb')); l = int(sys.argv[2]) - 1; ci0 = int(sys.argv[3]); kmax = int(sys.argv[4])
nal = D['nal']; P2 = alpha3.P2
T = D.get('T'); V = D.get('V')
if T is not None:
    ch = v4torus.root_choices(T)
    roots = list(itertools.product(*ch))[ci0]
    r_ = random.Random(5)
    al = v4torus.alpha_from(T, roots, [r_.randrange(1, P2) for _ in range(T['N'] - T['k'])])
    al = [x if x is not None else r_.randrange(1, P2) for x in al]
    parts = v4torus.alpha_parts(T, roots)
else:
    a, b = V['sols'][ci0]
    al = alpha3.alpha_point(V, nal, a, b, random.Random(5)); parts = None
cl, seedcols, cap = D['sys'][l]
rows, cols, ent, kl = cl
nr, nc, A = MS.build(cl, al, random.Random(7 + l), 'seed', seedcols, P2)
rowsof = {}; colsof = {}
for (r, c), v in A.items():
    if v:
        rowsof.setdefault(c, set()).add(r); colsof.setdefault(r, set()).add(c)
seed = next(iter(seedcols))
Rset = {0}; Cset = {seed}; frontier = {seed}
for k in range(1, kmax + 1):
    newR = set()
    for c in frontier:
        newR |= rowsof.get(c, set())
    newR -= Rset; Rset |= newR
    newC = set()
    for r in newR:
        newC |= colsof.get(r, set())
    newC -= Cset; Cset |= newC; frontier = newC
    sub = {(r, c): v for (r, c), v in A.items() if r in Rset}
    ncols = len(set(c for (r, c) in sub))
    ok, rk, st = sparsemod.consistent_sparse(dict(sub), {0: 1}, P2)
    print("k=%d: rows %d, columns %d -> %s (rank %s)" % (k, len(Rset), ncols, 'consistent' if ok else 'INCONSISTENT', rk), flush=True)
    if ok is False:
        # which alpha / free-parameter monomials occur on these rows
        inv_rows = {v: k_ for k_, v in rows.items()}
        ks = set()
        for (r, c, cf) in ent:
            if r + 1 in Rset:
                ks |= set(cf)
        print("   alphas on these rows: %s" % sorted(ks), flush=True)
        if parts is not None:
            print("   their free-parameter monomials: %s" % sorted(set(parts[k_][1] for k_ in ks if k_ in parts)), flush=True)
        pickle.dump(dict(Rset=Rset, k=k, l=l, ci=ci0), open('loc_%s_o%d_c%d.pkl' % (sys.argv[1].split('/')[-1][:-4], l + 1, ci0), 'wb'))
        break
