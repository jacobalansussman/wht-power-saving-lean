"""alpha3.py (w8-b5sparse): the d = 3 test ON the exact solution variety V_4 of the d = 4 system (no LM, exact mod p).

Why: the d = 4 passes of the orbit designs hold only at SPECIAL alpha (modsolve at a random alpha says INCONSISTENT at
d = 4 although LM reaches 1e-31).  Every d = 3 solution has its alpha in V_4 (the d = 3 system contains all d = 4 rows
with the same complete producer sets).  For the orbit designs the d = 4 closures have 2 live columns per output and
binomial rows  v_i alpha_i w_seed + v_j alpha_j w_other = 0.  With w_seed = 1 and s_l := -w_other(out l):
alpha_i = s_l (v_j / v_i) alpha_j.  Solving this graph gives  alpha_i = base(component) * C_i * s1^E1(i) * s2^E2(i)  and
cycle conditions s^v = 1 (v in Z^2); with all alpha != 0 the solutions (s1, s2) form a finite group (when the cycle
lattice has rank 2), enumerated exactly in F_p with p = 2305842973415982001 = 1 + k * lcm(1..24) (all roots of unity of
order <= 24 exist).  At d = 3 we then test exact consistency (FLINT rank) at every (s1, s2) solution and RANDOM bases:
the generic verdict on every component of V_4 with all alpha != 0.  (Components with some alpha = 0 are sub-designs.)
usage: as jd2.py (dlist 543)."""
import os, sys, time, random, itertools
from math import gcd
from fractions import Fraction
import jd2
import jdesign as J
import modsolve as MS

P2 = 2305842973415982001
ZL = 1465896768020624744
import os as _os
if _os.environ.get('PRIME2'):                       # second prime, same root-of-unity property (cross-check)
    P2 = 2305837608478644241; ZL = 558701126638716937
LL = 5354228880
V4 = {}
RELAX = int(os.environ.get('RELAX', '1'))


def zeta(N):
    assert LL % N == 0
    return pow(ZL, LL // N, P2)


def v4_binomial(stash, nal):
    edges = []
    for l in (0, 1):
        s = stash[l]
        rows, cols, ent, kl = s['cl']
        if len(cols) != 2 or len(s['seedcols']) != 1:
            return None, 'out %d: %d columns (need 2 with one seed)' % (l + 1, len(cols))
        seed = next(iter(s['seedcols'])); other = 1 - seed
        byrow = {}
        for (r, c, cf) in ent:
            byrow.setdefault(r, {}).setdefault(c, []).append(cf)
        for r, d in byrow.items():
            if set(d) != {0, 1} or any(len(v) != 1 or len(v[0]) != 1 for v in d.values()):
                return None, 'out %d: non-binomial row %s' % (l + 1, d)
            (ki, vi), = d[seed][0].items(); (kj, vj), = d[other][0].items()
            e = (1, 0) if l == 0 else (0, 1)
            edges.append((ki, kj, Fraction(vj).limit_denominator(10 ** 6) / Fraction(vi).limit_denominator(10 ** 6), e))
    # potentials: alpha_k = base * C_k * s^E_k ; node 0 is the constant 1
    adj = {}
    for (i, j, c, e) in edges:
        adj.setdefault(j, []).append((i, c, e))                 # alpha_i = c s^e alpha_j
        adj.setdefault(i, []).append((j, 1 / c, (-e[0], -e[1])))
    pot = {}; comp = {}; cyc = []
    for root in [0] + list(range(1, nal + 1)):
        if root in pot:
            continue
        pot[root] = (Fraction(1), (0, 0)); comp[root] = root
        stack = [root]
        while stack:
            j = stack.pop()
            Cj, Ej = pot[j]
            for (i, c, e) in adj.get(j, ()):
                Ci, Ei = c * Cj, (e[0] + Ej[0], e[1] + Ej[1])
                if i not in pot:
                    pot[i] = (Ci, Ei); comp[i] = root; stack.append(i)
                else:
                    Ci0, Ei0 = pot[i]
                    v = (Ei[0] - Ei0[0], Ei[1] - Ei0[1])
                    if Ci != Ci0:
                        if v == (0, 0):
                            return None, 'inconsistent constants (alpha forced 0)'
                        return None, 'constant cycle condition %s = %s' % (v, Ci0 / Ci)
                    if v != (0, 0):
                        cyc.append(v)
    # cycle lattice
    g1 = 0
    for a in cyc:
        for b in cyc:
            g1 = gcd(g1, abs(a[0] * b[1] - a[1] * b[0]))
    info = dict(ncyc=len(cyc), cycset=sorted(set(cyc))[:12], comps=sorted(set(comp.values())), order=g1)
    if g1 == 0:
        return None, 'cycle lattice rank < 2: %s' % info
    N = g1
    if LL % N:
        return None, 'group order %d does not divide lcm(1..24)' % N
    sols = [(a, b) for a in range(N) for b in range(N) if all((a * v[0] + b * v[1]) % N == 0 for v in cyc)]
    return dict(pot=pot, comp=comp, N=N, sols=sols, info=info), None


def alpha_point(V, nal, a, b, rng):
    z = zeta(V['N']); s1 = pow(z, a, P2); s2 = pow(z, b, P2)
    bases = {}
    al = []
    for k in range(1, nal + 1):
        if k not in V['pot']:
            al.append(rng.randrange(1, P2)); continue
        C, E = V['pot'][k]; rt = V['comp'][k]
        if rt not in bases:
            bases[rt] = 1 if rt == 0 else rng.randrange(1, P2)
        x = bases[rt] * (C.numerator % P2) * pow(C.denominator % P2, P2 - 2, P2) % P2
        x = x * pow(s1, E[0] % (P2 - 1), P2) % P2 * pow(s2, E[1] % (P2 - 1), P2) % P2
        al.append(x)
    return al


DENSEMAX = float(os.environ.get('DENSEMAX', '4e6'))


def consistent_at(cl, alpha, seed, norm, seedcols):
    import flint
    rng = random.Random(seed)
    nr, nc, A = MS.build(cl, alpha, rng, norm, seedcols, P2)
    if nc == 0:
        return False, 0, 0, nr, nc
    if nr * nc > DENSEMAX:
        import sparsemod
        ok, rk, st = sparsemod.consistent_sparse(A, {0: 1}, P2)
        if ok is None:
            print("        sparse elimination gave up (fill %s)" % st, flush=True)
            return True, rk, rk, nr, nc
        return ok, rk, rk + (0 if ok else 1), nr, nc
    Ma = flint.nmod_mat(nr, nc, P2); Mt = flint.nmod_mat(nr, nc + 1, P2)
    for (r, c), v in A.items():
        if v:
            Ma[r, c] = v; Mt[r, c] = v
    Mt[0, nc] = 1
    ra = Ma.rank(); rt = Mt.rank()
    return ra == rt, ra, rt, nr, nc


def relaxed_all_bases(cl, V, nal, a, b, seed, norm, seedcols):
    """EXACT test for ALL bases at once on the V4 component (a, b): alpha_k = base(comp k) * const_k, so
    N(alpha) = N_0 + sum_comp base_comp * N_comp.  Writing u_comp = base_comp * w, a solution gives
    N_0 w + sum N_comp u_comp = e: so if e is NOT in colspan [N_0 | N_comp ...], no choice of bases works."""
    import flint
    rows, cols, ent, kl = cl
    klass, kterms = kl
    z = zeta(V['N']); s1 = pow(z, a, P2); s2 = pow(z, b, P2)
    const = {}; part = {}
    for k in range(1, nal + 1):
        if k not in V['pot']:
            part[k] = ('free', k); const[k] = 1; continue
        C, E = V['pot'][k]; rt = V['comp'][k]
        x = (C.numerator % P2) * pow(C.denominator % P2, P2 - 2, P2) % P2
        x = x * pow(s1, E[0] % (P2 - 1), P2) % P2 * pow(s2, E[1] % (P2 - 1), P2) % P2
        const[k] = x; part[k] = rt if rt != 0 else 0
    parts = sorted(set(part.values()) | {0}, key=str)
    pidx = {pp: i for i, pp in enumerate(parts)}
    nr, nc = len(rows) + 1, len(cols)
    A = {}
    def add(r, c, cf, mult=1):
        for k, v in cf.items():
            pp = 0 if k == 0 else part[k]
            x = MS.fval(v, P2) * (1 if k == 0 else const[k]) % P2 * mult % P2
            key = (r, pidx[pp] * nc + c); A[key] = (A.get(key, 0) + x) % P2
    for (r, c, cf) in ent:
        add(r + 1, c, cf)
    rng = random.Random(seed)
    seedrows = None
    if norm == 'seed':
        seedrows = set(rk for c in seedcols for (rk, cf) in kterms.get(c, []))
    lam = {}
    for c in klass:
        for (rk, cf) in kterms.get(c, []):
            if seedrows is not None and rk not in seedrows:
                continue
            if rk not in lam:
                lam[rk] = rng.randrange(1, P2)
            add(0, c, cf, lam[rk])
    W = len(parts) * nc
    if nr * W > DENSEMAX:
        import sparsemod
        ok, rk, st = sparsemod.consistent_sparse(A, {0: 1}, P2)
        if ok is None:
            return True, rk, rk, nr, W, len(parts)
        return ok, rk, rk + (0 if ok else 1), nr, W, len(parts)
    Ma = flint.nmod_mat(nr, W, P2); Mt = flint.nmod_mat(nr, W + 1, P2)
    for (r, c), v in A.items():
        if v:
            Ma[r, c] = v; Mt[r, c] = v
    Mt[0, W] = 1
    ra = Ma.rank(); rt = Mt.rank()
    return ra == rt, ra, rt, nr, W, len(parts)


def relaxed_parts(cl, nal, partof, const, seed, seedcols, norm='seed'):
    """exact relaxed test with alpha_k = const_k * monomial_{partof[k]}: e in colspan [N_part ...] ?"""
    import flint
    rows, cols, ent, kl = cl
    klass, kterms = kl
    parts = sorted(set(partof.values()) | {()}, key=str)
    pidx = {pp: i for i, pp in enumerate(parts)}
    nr, nc = len(rows) + 1, len(cols)
    A = {}
    def add(r, c, cf, mult=1):
        for k, v in cf.items():
            pp = () if k == 0 else partof[k]
            x = MS.fval(v, P2) * (1 if k == 0 else const[k]) % P2 * mult % P2
            key = (r, pidx[pp] * nc + c); A[key] = (A.get(key, 0) + x) % P2
    for (r, c, cf) in ent:
        add(r + 1, c, cf)
    rng = random.Random(seed)
    seedrows = set(rk for c in seedcols for (rk, cf) in kterms.get(c, []))
    lam = {}
    for c in klass:
        for (rk, cf) in kterms.get(c, []):
            if norm == 'seed' and rk not in seedrows:
                continue
            if rk not in lam:
                lam[rk] = rng.randrange(1, P2)
            add(0, c, cf, lam[rk])
    W = len(parts) * nc
    import sparsemod
    if nr * W > DENSEMAX:
        ok, rk, st = sparsemod.consistent_sparse(A, {0: 1}, P2)
        if ok is None:
            return True, rk, rk, nr, W, len(parts)
        return ok, rk, rk + (0 if ok else 1), nr, W, len(parts)
    Ma = flint.nmod_mat(nr, W, P2); Mt = flint.nmod_mat(nr, W + 1, P2)
    for (r, c), v in A.items():
        if v:
            Ma[r, c] = v; Mt[r, c] = v
    Mt[0, W] = 1
    ra = Ma.rank(); rt = Mt.rank()
    return ra == rt, ra, rt, nr, W, len(parts)


def torus_d3(T, nal, d):
    import v4torus
    ch = v4torus.root_choices(T)
    if jd2.STASH[0]['stats'].get('cap') or jd2.STASH[1]['stats'].get('cap'):
        import pickle
        pickle.dump(dict(T=T, nal=nal, d=d, sys=[(jd2.STASH[l]['cl'], jd2.STASH[l]['seedcols'], jd2.STASH[l]['stats'].get('cap')) for l in (0, 1)]),
                    open('sys_d%d_%s.pkl' % (d, sys.argv[1]), 'wb'))
        if os.environ.get('PICKLEONLY'):
            print('     d=%d TRUNCATED system pickled (PICKLEONLY)' % d, flush=True); return 1.0, None
        nex = 0; nall = 0
        for ci, roots in enumerate(itertools.product(*ch)):
            pa = v4torus.alpha_parts(T, roots)
            partof = {}; const = {}
            for k in range(1, nal + 1):
                if k in pa:
                    const[k] = pa[k][0]; partof[k] = pa[k][1]
                else:
                    const[k] = 1; partof[k] = ('free', k)
            rr = []
            for l in (0, 1):
                s = jd2.STASH[l]
                t = time.time()
                x = relaxed_parts(s['cl'], nal, partof, const, 9 + l, s['seedcols'], 'seed')
                rr.append(x[0])
                print("     d=%d TRUNCATED, torus root choice %d, relaxed out %d (seed norm): %d x %d (%d parts), rank %d/%d -> %s (%.0fs)" % (
                    d, ci, l + 1, x[3], x[4], x[5], x[1], x[2], 'consistent (inconclusive)' if x[0] else 'EXCLUDED', time.time() - t), flush=True)
            nex += (not all(rr)); nall += 1
        print("     d=%d verdict on V4 (torus) [truncated system]: the relaxed test excludes %d of %d root choices (a 'pass' below means INCONCLUSIVE)" % (d, nex, nall), flush=True)
        return (0.0 if nex < nall else 1.0), None
    anyok = False; nall = 0; nex = 0
    for ci, roots in enumerate(itertools.product(*ch)):
        if ci >= 256:
            print("     (torus: more than 256 root choices; stopped at 256 -- NOT exhaustive)", flush=True)
            break
        r_ = random.Random(1000 + ci)
        al = v4torus.alpha_from(T, roots, [r_.randrange(1, P2) for _ in range(T['N'] - T['k'])])
        al = [x if x is not None else r_.randrange(1, P2) for x in al]
        okl = []
        for l in (0, 1):
            s = jd2.STASH[l]
            okl.append(consistent_at(s['cl'], al, 7 + l, 'seed', s['seedcols']))
        both = all(x[0] for x in okl)
        anyok |= both; nall += 1
        line = "     d=%d on V4 (torus) root choice %d: %s" % (d, ci, ' '.join('o%d:%s(%d/%d)' % (l + 1, 'C' if x[0] else 'I', x[1], x[2]) for l, x in enumerate(okl)))
        if not both and RELAX:
            pa = v4torus.alpha_parts(T, roots)
            partof = {}; const = {}
            for k in range(1, nal + 1):
                if k in pa:
                    const[k] = pa[k][0]; partof[k] = pa[k][1]
                else:
                    const[k] = 1; partof[k] = ('free', k)
            rr = []
            for l in (0, 1):
                s = jd2.STASH[l]
                rr.append(relaxed_parts(s['cl'], nal, partof, const, 9 + l, s['seedcols'], 'all'))
            ex = any(not x[0] for x in rr)
            nex += ex
            line += "; relaxed all-parameters: %s" % ' '.join('o%d:%s(%d/%d, %d parts)' % (l + 1, 'EXCL' if not x[0] else 'incl', x[1], x[2], x[5]) for l, x in enumerate(rr))
        print(line, flush=True)
    st = jd2.STASH[0]['stats']
    print("     d=%d verdict on V4 (torus; all unknowns of the seed component nonzero): %s; relaxed test excludes %d of %d root choices%s" % (
        d, 'some point CONSISTENT' if anyok else 'INCONSISTENT at every root choice (generic free parameters)', nex, nall,
        ' [truncated system]' if st.get('cap') else ' [complete closure]'), flush=True)
    return (0.0 if anyok else 1.0), None


def lm3(Nss, nal, rng, starts=3, iters=60):
    d = jd2.STASH[0]['d']
    if d >= 5:
        return jd2.lm2(Nss, nal, rng, starts, iters)
    if d == 4:
        V, why = v4_binomial(jd2.STASH, nal)
        if V is None:
            import v4torus
            T, why2 = v4torus.solve(jd2.STASH, nal)
            if T is None:
                print("     V4: not handled (%s; torus: %s); falling back to jd2" % (why, why2), flush=True)
                return jd2.lm2(Nss, nal, rng, starts, iters)
            V4['T'] = T
            ch = v4torus.root_choices(T)
            ncomb = 1
            for c_ in ch:
                ncomb *= len(c_)
            print("     V4 (torus, exact): %d unknowns, %d equations, rank %d, elementary divisors %s, sigma %s, %d free "
                  "parameters, %d root choices" % (T['N'], T['neq'], T['k'], T['ds'], [str(x) for x in T['sig']],
                                                    T['N'] - T['k'], ncomb), flush=True)
            bad = 0; nchk = 0
            for roots in itertools.islice(itertools.product(*ch), 16):
                r_ = random.Random(nchk); nchk += 1
                al = v4torus.alpha_from(T, roots, [r_.randrange(1, P2) for _ in range(T['N'] - T['k'])])
                al = [x if x is not None else r_.randrange(1, P2) for x in al]
                for l in (0, 1):
                    s = jd2.STASH[l]
                    bad += (not consistent_at(s['cl'], al, 3 + l, 'seed', s['seedcols'])[0])
            print("     V4 (torus) check: %d root choices consistent at d=4 exactly: %s" % (nchk, bad == 0), flush=True)
            return (0.0 if bad == 0 else 1.0), None
        V4['V'] = V
        print("     V4 (exact): alpha_k = base * C_k * s1^E1 s2^E2; components %s; cycle lattice order %d; %d solutions "
              "(s1, s2) = (z^a, z^b), z = zeta_%d: %s" % (V['info']['comps'], V['N'], len(V['sols']), V['N'], V['sols']), flush=True)
        # check: each solution makes d=4 consistent exactly
        bad = 0
        for (a, b) in V['sols']:
            al = alpha_point(V, nal, a, b, random.Random(5))
            for l in (0, 1):
                s = jd2.STASH[l]
                ok = consistent_at(s['cl'], al, 3 + l, 'seed', s['seedcols'])[0]
                bad += (not ok)
        print("     V4 check: every solution consistent at d=4 exactly: %s" % (bad == 0), flush=True)
        return (0.0 if bad == 0 else 1.0), None
    V = V4.get('V')
    if V is None and 'T' in V4:
        return torus_d3(V4['T'], nal, d)
    if V is None:
        ncmax = max(Ns[0].shape[1] for Ns in Nss)
        if ncmax <= jd2.CFG['LMMAX']:
            return jd2.lm2(Nss, nal, rng, starts, iters)
        # big system, V4 unknown: only the GENERIC-alpha exact test is available (special alpha NOT covered)
        oks = []
        for l in (0, 1):
            s = jd2.STASH[l]
            r_ = random.Random(77 + l)
            al = [r_.randrange(1, P2) for _ in range(nal)]
            ok, ra, rt, nr, nc = consistent_at(s['cl'], al, 7 + l, 'seed', s['seedcols'])
            oks.append(ok)
            print("     d=%d GENERIC-ALPHA ONLY (V4 not binomial) out %d: %d x %d, rank %d/%d -> %s%s" % (
                d, l + 1, nr, nc, ra, rt, 'consistent' if ok else 'inconsistent',
                ' [truncated system]' if s['stats'].get('cap') else ''), flush=True)
        return (0.0 if all(oks) else 1.0), None
    anyok = False
    if jd2.STASH[0]['stats'].get('cap') or jd2.STASH[1]['stats'].get('cap'):
        import pickle
        pickle.dump(dict(V=V, nal=nal, d=d, sys=[(jd2.STASH[l]['cl'], jd2.STASH[l]['seedcols'], jd2.STASH[l]['stats'].get('cap')) for l in (0, 1)]),
                    open('sys_d%d_%s.pkl' % (d, sys.argv[1]), 'wb'))
        if os.environ.get('PICKLEONLY'):
            print('     d=%d TRUNCATED system pickled (PICKLEONLY)' % d, flush=True); return 1.0, None
        nex = 0
        for (a, b) in V['sols']:
            rr = []
            for l in (0, 1):
                s = jd2.STASH[l]
                t = time.time()
                ok, ra, rt, nr, W, npart = relaxed_all_bases(s['cl'], V, nal, a, b, 9 + l, 'seed', s['seedcols'])
                rr.append(ok)
                print("     d=%d TRUNCATED, V4 point (a,b)=(%d,%d), relaxed all-bases out %d (seed norm): %d x %d, rank %d/%d -> %s (%.0fs)" % (
                    d, a, b, l + 1, nr, W, ra, rt, 'consistent (inconclusive)' if ok else 'EXCLUDED', time.time() - t), flush=True)
            nex += (not all(rr))
        print("     d=%d verdict on V4 [truncated system]: the relaxed test excludes %d of %d root choices (a 'pass' below means INCONCLUSIVE)" % (d, nex, len(V['sols'])), flush=True)
        return (0.0 if nex < len(V['sols']) else 1.0), None
    for (a, b) in V['sols']:
        res = []
        for trial in range(2):
            al = alpha_point(V, nal, a, b, random.Random(100 * a + 10 * b + trial))
            for l in (0, 1):
                s = jd2.STASH[l]
                for norm in ('seed', 'all'):
                    t = time.time()
                    ok, ra, rt, nr, nc = consistent_at(s['cl'], al, 7 + l, norm, s['seedcols'])
                    res.append((l + 1, norm, ok, ra, rt))
        both = [all(x[2] for x in res if x[1] == norm and True) for norm in ('seed', 'all')]
        print("     d=%d on V4 point (a,b)=(%d,%d): %s -> %s" % (
            d, a, b, ' '.join('o%d/%s:%s(%d/%d)' % (o, n, 'C' if ok else 'I', ra, rt) for (o, n, ok, ra, rt) in res[:4]),
            'CONSISTENT (seed norm, both outputs)' if both[0] else 'inconsistent'), flush=True)
        anyok |= both[0]
        if not both[0] and RELAX:
            rr = []
            for l in (0, 1):
                s = jd2.STASH[l]
                t = time.time()
                for nm in ('seed', 'all'):
                    ok, ra, rt, nr, W, npart = relaxed_all_bases(s['cl'], V, nal, a, b, 9 + l, nm, s['seedcols'])
                    rr.append(ok)
                    print("        relaxed ALL-bases test out %d (%s norm): %d x %d (%d parts), rank %d / %d -> %s (%.0fs)" % (
                        l + 1, nm, nr, W, npart, ra, rt, 'inconclusive (consistent)' if ok else 'EXCLUDED for every choice of bases', time.time() - t), flush=True)
    st = jd2.STASH[0]['stats']
    print("     d=%d verdict on V4 (all alpha != 0, generic bases): %s%s" % (
        d, 'some point CONSISTENT' if anyok else 'INCONSISTENT at every solution', ' [truncated system]' if st.get('cap') else ' [complete closure]'), flush=True)
    return (0.0 if anyok else 1.0), None


J.lm = lm3

if __name__ == "__main__":
    print("alpha3 config %s argv %s PYTHONHASHSEED=%s" % (jd2.CFG, sys.argv[1:], os.environ.get('PYTHONHASHSEED')), flush=True)
    J.main()
