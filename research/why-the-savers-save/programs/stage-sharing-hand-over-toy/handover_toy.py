"""s-wide / handover_toy.py -- an exact two-stage hand-over toy in the WIDER game (additions need no common frame).

Model (exact engine semantics, restricted to the sub-game of directions + shifts; no phases, no address maps):
  * Block size n = 3.  A frame is q : F2^3 -> Z4 (Fourier side: dir(z) multiplies character y by i^<y,z>), so
    dir(z) adds chi_z(y) = <y,z> to q, shift(w) adds 2*chi_w (free).  Frames are compared modulo 2*linear (shifts).
  * Bookkeeping: input l starts at frame s_l (x_T at chi_{a_T}, the line of port T; y_T at 0).  The group is abelian,
    so bookkeeping is a translation and moves act on every term of a register alike.  A register may hold terms at
    several frames (wider game); terms merge or cancel only when (input, frame) are identical.
  * Toy ports: three mutually orthogonal labels a = e1, e2, e3 in U = F2^3 (each knows the other two).
    Stage i (fwd):   X_T climbs chi_a -> chi_a + chi_d1 -> F = chi_e1 + chi_e2 + chi_e3 ; Y_T climbs 0 -> chi_d1 -> psi_T
                     (psi_T = chi_d1 + chi_d2, D_T = the other two unit vectors); stage-i scatter: Y_T += x_T at psi_T,
                     done in stage i by a copy of x_T moved F -> psi_T (1 move; paid in stage i).
    Stage i+1 (rev): X_S -= Y_S, i.e. X_S at frame F must receive (y_S + x_S) at F.
  * A design for stage i+1 is a forest of scratch registers: every edge is one paid move (a direction z), content may be
    added at any node (free, any frame), a copy at a node duplicates the register into its children, and every leaf is
    added into the target.  The sum over leaves must equal the target exactly (no junk left in a live array).
    Feasibility for a given forest = a linear system over Q (solved mod p = 2^61-1; generic, see report).
  * Sources N (no hand-over): stage-(i+1) data only: x_T at F (X_T), (y_T + x_T) at psi_T (Y_T after stage i).
    Sources H (hand-over): N plus everything stage i had for free or already paid: y_T at 0, chi_d1, psi_T (Y_T before
    the stage-i scatter), x_T at chi_a, chi_a + chi_d1 (X_T's path), and the stage-i copies x_T at psi_T.
  Output: the least number of paid moves for each target and source set, by exhaustive enumeration of forests with up to
  K edges (K = 4 by default).   usage: python3 -I handover_toy.py [K]
"""
import itertools, sys
sys.dont_write_bytecode = True
N = 3; Y = range(1 << N); P = (1 << 61) - 1

def chi(z): return tuple(bin(y & z).count('1') & 1 for y in Y)
def add(q, r): return tuple((a + b) & 3 for a, b in zip(q, r))
def canon(q):  # modulo free shifts (2*chi_w)
    return min(add(q, tuple(2 * c for c in chi(w))) for w in Y)
ZERO = canon(tuple(0 for _ in Y))
DIRS = list(range(1, 1 << N))
E = [1, 2, 4]
F = canon(add(add(chi(1), chi(2)), chi(4)))
ports = {a: [d for d in E if d != a] for a in E}          # D_T in climb order
psi = {a: canon(add(chi(ports[a][0]), chi(ports[a][1]))) for a in E}
d1f = {a: canon(chi(ports[a][0])) for a in E}

def src_N():
    S = []
    for a in E:
        S.append(({('x', a): 1}, F))
        S.append(({('y', a): 1, ('x', a): 1}, psi[a]))
    return S
def src_H():
    S = src_N()
    for a in E:
        S += [({('y', a): 1}, ZERO), ({('y', a): 1}, d1f[a]), ({('y', a): 1}, psi[a])]
        S += [({('x', a): 1}, canon(chi(a))), ({('x', a): 1}, canon(add(chi(a), chi(ports[a][0])))),
              ({('x', a): 1}, psi[a])]
    return S

def rank_mod(rows, ncols):
    rows = [r[:] for r in rows]; rk = 0
    for c in range(ncols):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % P), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], P - 2, P)
        rows[rk] = [v * inv % P for v in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % P:
                f = rows[i][c]; rows[i] = [(a - f * b) % P for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def feasible(nodes_disp, sources, target):
    """nodes_disp: for each node, list of displacements to the leaves below it (a copy duplicates).  Unknown
    lambda_{s,v}; equation per (input, final frame): sum lambda * coeff = target."""
    keys = {}; cols = []
    for v, disps in enumerate(nodes_disp):
        for (c, f) in sources:
            col = {}
            for dsp in disps:
                g = canon(add(f, dsp))
                for inp, co in c.items():
                    k = (inp, g); col[k] = col.get(k, 0) + co
            cols.append(col)
    allk = set(target)
    for col in cols: allk |= set(col)
    allk = sorted(allk); idx = {k: i for i, k in enumerate(allk)}
    A = [[0] * (len(cols) + 1) for _ in allk]
    for j, col in enumerate(cols):
        for k, co in col.items(): A[idx[k]][j] = co % P
    for k, co in target.items(): A[idx[k]][-1] = co % P
    rA = rank_mod([r[:-1] for r in A], len(cols)); rAb = rank_mod(A, len(cols) + 1)
    return rA == rAb

def forests(K):
    """Rooted forests with K edges, given as parent arrays over non-root nodes; roots are nodes with parent None.
    We enumerate trees edge by edge: node 0..; every new node is a child of an earlier node or a new root."""
    def rec(par, edges):
        if edges == K: yield list(par); return
        n = len(par)
        # new node as child of an existing node (one edge)
        for p in range(n):
            par.append(p); yield from rec(par, edges + 1); par.pop()
        # a new root (no edge), only if the last node is not a childless root (avoid duplicates of empty roots)
        if n == 0 or not (par[-1] is None):
            par.append(None); yield from rec(par, edges); par.pop()
    yield from rec([None], 0)

def designs(K):
    seen = set()
    for par in forests(K):
        ch = {i: [j for j, p in enumerate(par) if p == i] for i in range(len(par))}
        leaves = [i for i in range(len(par)) if not ch[i]]
        key = tuple(par)
        if key in seen: continue
        seen.add(key)
        edges = [i for i, p in enumerate(par) if p is not None]
        for labs in itertools.product(DIRS, repeat=len(edges)):
            lab = dict(zip(edges, labs))
            def disp_to(v, leaf):  # displacement from node v down to leaf (v ancestor of leaf)
                d = ZERO; u = leaf
                while u != v:
                    if par[u] is None: return None
                    d = canon(add(d, chi(lab[u]))); u = par[u]
                return d
            nodes_disp = []
            for v in range(len(par)):
                ds = [disp_to(v, l) for l in leaves]
                nodes_disp.append([d for d in ds if d is not None])
            yield len(edges), nodes_disp, (par, lab)

def least(target, sources, K):
    if feasible([[ZERO]], sources, target): return 0, None
    for k in range(1, K + 1):
        for ne, nd, dsg in designs(k):
            if ne != k: continue
            if feasible(nd, sources, target): return k, dsg
    return None, None

def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    print('frames: F =', F, '; psi =', psi, '; distinct psi:', len(set(psi.values())))
    T1 = {(('y', 1), F): 1, (('x', 1), F): 1}
    T2 = {}
    for a in E:
        T2[(('y', a), F)] = 1; T2[(('x', a), F)] = 1
    T2y = {(('y', a), F): 1 for a in E}
    for name, tg in [('T1 one port: (y1+x1) at F', T1), ('T2 coded total: sum_T (yT+xT) at F', T2),
                     ('T2y the y-part only: sum_T yT at F', T2y)]:
        for sname, S in [('N (no hand-over)', src_N()), ('H (hand-over, wider)', src_H())]:
            k, d = least(tg, S, K)
            print('%-38s | %-22s | least paid moves: %s%s' % (name, sname, k if k is not None else '> %d' % K,
                  '' if d is None else '   design parents=%s dirs=%s' % (d[0], d[1])))
            sys.stdout.flush()

if __name__ == '__main__':
    main()
