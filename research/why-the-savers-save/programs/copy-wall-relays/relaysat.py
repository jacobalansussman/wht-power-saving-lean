"""relaysat.py (d-relay) -- exact relay minrank over F_2 by SAT (python-sat Cadical153).

Model.  B = X Y^T over F_2 with X, Y in F_2^{v x r}.  B[S,S] = 1.  Ports carry a layer in {0..L-1}.  For S != T:
  - S knows T                      : B[S,T] free (S cancels x_T at a shared frame);
  - S does not know T              : B[S,T] = 0 unless layer(T) < layer(S)  (S relays T: reads T's register y_T + x_T).
Closure modes (used(S,T) := B[S,T] = 1 with S not knowing T):
  none  : decoder-level relaxation only (g-recombine/relay.py's model);
  one   : used(S,T) -> B[T,T'] = 0 for every T' known to T but not to S   (the brief's clause; T' = S included);
  trans : footprint F(T) = {T' known to T : B[T,T'] = 1} u U_{used(T,T'')} F(T'');
          used(S,T) -> F(T) inside known(S)  (labels cancelled in T's past must lie in the relay frame G <= a_S^perp).
Gauge (optional, --gauge): an 8-clique K (pairwise not knowing) forced into ONE common layer; then B[K,K] = I and
X_K = Y_K = I without loss of generality inside that restriction.
Gauge --gaugeX: X_K = I only, K unrestricted.  WITHOUT LOSS OF GENERALITY for every layering and every closure mode:
B[K,K] is unitriangular in the layer order (B[k,k'] = 0 unless layer(k') < layer(k)), hence invertible, so a rank-r
(r = |K|) B = X Y^T has X_K invertible and (X, Y) -> (X X_K^-1, Y X_K^T) keeps B.
usage: venv/bin/python -I relaysat.py INSTANCE r L MODE [--gauge] [--timeout s] [--save file.json]"""
import sys, time, json, itertools, threading
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from patterns import instance, max_clique_nonknow
from pysat.solvers import Cadical153, Solver
from pysat.formula import IDPool

def build(knows, r, L, mode, gauge=None, xonly=False):
    v = len(knows)
    pool = IDPool(); cl = []
    T_, F_ = 'T', 'F'
    X = [[pool.id(('x', S, k)) for k in range(r)] for S in range(v)]
    Y = [[pool.id(('y', S, k)) for k in range(r)] for S in range(v)]
    if gauge:
        for i, kk in enumerate(gauge):
            for k in range(r):
                X[kk][k] = T_ if k == i else F_
                if not xonly: Y[kk][k] = T_ if k == i else F_
    def AND(a, b):
        if a == F_ or b == F_: return F_
        if a == T_: return b
        if b == T_: return a
        z = pool.id(('and', a, b)); cl.extend([[-z, a], [-z, b], [z, -a, -b]]); return z
    def XOR(lits):
        const = 0; ls = []
        for l in lits:
            if l == T_: const ^= 1
            elif l == F_: pass
            else: ls.append(l)
        if not ls: return T_ if const else F_
        acc = ls[0]
        for l in ls[1:]:
            z = pool.id(('xor', acc, l)); cl.extend([[-z, acc, l], [-z, -acc, -l], [z, -acc, l], [z, acc, -l]]); acc = z
        if const:
            z = pool.id(('not', acc)); cl.extend([[z, acc], [-z, -acc]]); acc = z
        return acc
    Bc = {}
    def Bv(S, T):
        if (S, T) not in Bc: Bc[S, T] = XOR([AND(X[S][k], Y[T][k]) for k in range(r)])
        return Bc[S, T]
    def lit_true(l):   # clause meaning l holds
        if l == T_: return None
        if l == F_: return []
        return [l]
    def neg(l):
        return F_ if l == T_ else T_ if l == F_ else -l
    def clause(lits):
        out = []
        for l in lits:
            if l == T_: return
            if l == F_: continue
            out.append(l)
        cl.append(out)
    # layers: order encoding ge[S][l] = layer(S) >= l, l = 1..L-1
    ge = [[pool.id(('ge', S, l)) for l in range(1, L)] for S in range(v)]
    for S in range(v):
        for l in range(len(ge[S]) - 1): cl.append([-ge[S][l + 1], ge[S][l]])
    if gauge and L > 1 and not xonly:
        k0 = gauge[0]
        for kk in gauge[1:]:
            for l in range(L - 1): cl.extend([[-ge[kk][l], ge[k0][l]], [ge[kk][l], -ge[k0][l]]])
    def before(T, S):
        """literal list: clause part meaning layer(T) < layer(S) (returns a fresh var p with p -> T before S)"""
        if L == 1: return F_
        p = pool.id(('bef', T, S)); zs = []
        for l in range(L - 1):
            z = pool.id(('z', T, S, l)); cl.extend([[-z, -ge[T][l]], [-z, ge[S][l]]]); zs.append(z)
        cl.append([-p] + zs)
        return p
    for S in range(v):
        clause([Bv(S, S)])
        for T in range(v):
            if T == S or knows[S][T]: continue
            clause([neg(Bv(S, T)), before(T, S)])
    nonknow = [(S, T) for S in range(v) for T in range(v) if S != T and not knows[S][T]]
    if mode == 'one':
        for S, T in nonknow:
            for Tp in range(v):
                if Tp != T and knows[T][Tp] and not knows[S][Tp]:
                    clause([neg(Bv(S, T)), neg(Bv(T, Tp))])
    elif mode == 'trans':
        f = {(T, Tp): pool.id(('f', T, Tp)) for T in range(v) for Tp in range(v) if Tp != T and knows[T][Tp]}
        for (T, Tp), fv in f.items(): clause([neg(Bv(T, Tp)), fv])
        for T, Tpp in nonknow:            # T relays Tpp: F(Tpp) must sit inside known(T), and joins F(T)
            for Tp in range(v):
                if (Tpp, Tp) not in f: continue
                if Tp != T and knows[T][Tp]: clause([neg(Bv(T, Tpp)), -f[Tpp, Tp], f[T, Tp]])
                else: clause([neg(Bv(T, Tpp)), -f[Tpp, Tp]])
    elif mode != 'none':
        raise SystemExit('mode?')
    return pool, cl, X, Y, ge

def decode(model, X, Y, ge, v, r, L):
    m = set(l for l in model if l > 0)
    val = lambda l: True if l == 'T' else False if l == 'F' else (l in m)
    Xv = [[int(val(X[S][k])) for k in range(r)] for S in range(v)]
    Yv = [[int(val(Y[S][k])) for k in range(r)] for S in range(v)]
    lay = [sum(1 for l in ge[S] if l in m) for S in range(v)]
    return Xv, Yv, lay

def verify(knows, Xv, Yv, lay, mode):
    """independent check; returns dict of findings"""
    v = len(knows); r = len(Xv[0])
    B = [[sum(Xv[S][k] & Yv[T][k] for k in range(r)) & 1 for T in range(v)] for S in range(v)]
    bad_diag = sum(1 for S in range(v) if B[S][S] != 1)
    used = [(S, T) for S in range(v) for T in range(v) if S != T and not knows[S][T] and B[S][T]]
    bad_order = sum(1 for S, T in used if not lay[T] < lay[S])
    viol_one = sum(1 for S, T in used for Tp in range(v) if Tp != T and knows[T][Tp] and not knows[S][Tp] and B[T][Tp])
    # transitive footprints in layer order
    Fp = {}
    for T in sorted(range(v), key=lambda t: lay[t]):
        s = {Tp for Tp in range(v) if Tp != T and knows[T][Tp] and B[T][Tp]}
        for (S2, T2) in used:
            if S2 == T: s |= Fp[T2]
        Fp[T] = s
    viol_trans = sum(1 for S, T in used for Tp in Fp[T] if not knows[S][Tp])
    # static re-derivation c^T (the no-gain construction): c^T = B[T] - sum_{used(T,T'')} B[T,T''] c^{T''}
    C = {}
    for T in sorted(range(v), key=lambda t: lay[t]):
        c = list(B[T])
        for (S2, T2) in used:
            if S2 == T: c = [(a ^ b) for a, b in zip(c, C[T2])]
        C[T] = c
    static_ok = all(C[S][S] == 1 and all(C[S][T] == 0 for T in range(v) if T != S and not knows[S][T]) for S in range(v))
    rkC = rank2([sum(b << i for i, b in enumerate(C[S])) for S in range(v)])
    rkB = rank2([sum(b << i for i, b in enumerate(B[S])) for S in range(v)])
    return dict(bad_diag=bad_diag, bad_order=bad_order, relays_used=len(used), viol_one=viol_one, viol_trans=viol_trans,
                rankB=rkB, cT_static_fit=static_ok, rank_cT=rkC, layer_sizes=[lay.count(l) for l in range(max(lay) + 1)])

def rank2(rows):
    basis = []
    for x in rows:
        for b in basis: x = min(x, x ^ b)
        if x: basis.append(x); basis.sort(reverse=True)
    return len(basis)

SOLVER = 'cd153'
def solve(name, r, L, mode, gauge=False, timeout=None, save=None, quiet=False, xonly=False):
    knows, G = instance(name); v = len(knows)
    K = None
    if gauge:
        K = max_clique_nonknow(knows, size=r)
        assert len(K) <= r
        if len(K) < r: assert xonly, 'partial clique needs --gaugeX (X_K = [I 0])'
    t0 = time.time()
    pool, cl, X, Y, ge = build(knows, r, L, mode, K, xonly)
    s = Cadical153(bootstrap_with=cl) if SOLVER == 'cd153' else Solver(name=SOLVER, bootstrap_with=cl)
    timer = None
    if timeout:
        timer = threading.Timer(timeout, lambda: s.interrupt()); timer.start()
        res = s.solve_limited(expect_interrupt=True)
    else:
        res = s.solve()
    if timer: timer.cancel()
    dt = time.time() - t0
    tag = '%s r=%d L=%d closure=%s gauge=%s' % (name, r, L, mode, ('XK=I(|K|=%d)' % len(K) if xonly else 'XK=YK=I,K-one-layer') if K else 'none')
    if res is None:
        print('%-60s UNDECIDED (timeout %ss)  vars %d clauses %d  %.1fs' % (tag, timeout, pool.top, len(cl), dt), flush=True)
        return None, None
    if not res:
        print('%-60s UNSAT  vars %d clauses %d  %.1fs' % (tag, pool.top, len(cl), dt), flush=True)
        return False, None
    Xv, Yv, lay = decode(s.get_model(), X, Y, ge, v, r, L)
    chk = verify(knows, Xv, Yv, lay, mode)
    print('%-60s SAT    vars %d clauses %d  %.1fs  check %s' % (tag, pool.top, len(cl), dt, chk), flush=True)
    if save:
        json.dump(dict(instance=name, r=r, L=L, mode=mode, gauge=K, X=Xv, Y=Yv, layer=lay, check=chk), open(save, 'w'))
    return True, chk

if __name__ == '__main__':
    a = sys.argv[1:]
    name, r, L, mode = a[0], int(a[1]), int(a[2]), a[3]
    gauge = '--gauge' in a or '--gaugeX' in a
    xonly = '--gaugeX' in a
    SOLVER = a[a.index('--solver') + 1] if '--solver' in a else 'cd153'
    timeout = float(a[a.index('--timeout') + 1]) if '--timeout' in a else None
    save = a[a.index('--save') + 1] if '--save' in a else None
    solve(name, r, L, mode, gauge, timeout, save, xonly=xonly)
