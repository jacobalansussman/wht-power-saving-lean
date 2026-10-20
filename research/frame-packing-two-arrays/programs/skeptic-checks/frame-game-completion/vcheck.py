"""vcheck.py (v-pk2-complete): independent check of pk2-complete Theorem 1 / Corollary 2 on pk2's programs.
usage: vcheck.py m fmax [modes]
 1. frame-game identity of each exported program, re-verified with MY Sp code over Q (my T_z = gbar(dir z),
    my Kbar = gbar(K), their p's converted from row bitmasks and checked to lie in P)
 2. lift with MY Lambda (veng.D; P lifts from my table of legitimate free-letter words; final correction
    Q = D(K^-1), which is what Lambda(Kbar^-1) gives -- a different formula from pk2's [K] D(K)^-1 [K]^-1)
 3. FORMAL check in C[G] with exact Gaussian-rational coefficients: output for input l == [K] (one term, coef 1)
 4. raw complex arrays at f = 1..fmax (every group element applied as g^(x)f)
 5. legitimacy / cost: every free group element = (legit word in P table) * i^j P_u; paid letters = the program's moves
 6. negative control: fp-build's f = 1 lift alone (mode f1)."""
import sys, json, time
import numpy as np
from fractions import Fraction as Fr
from veng import *

m = int(sys.argv[1]); fmax = int(sys.argv[2]); modes = sys.argv[3].split(',') if len(sys.argv) > 3 else ['full', 'f1']
E = VEng(m); E.build_Ptable(); n = E.n
progs = json.load(open(sys.argv[4] if len(sys.argv) > 4 else f'progs_m{m}.json'))

def conv(rows):
    M = np.zeros((2 * m, 2 * m), dtype=np.int64)
    for i, r in enumerate(rows):
        for j in range(2 * m): M[i, j] = (r >> j) & 1
    return M
Tz = {z: E.gbar(E.dirm(z)) for z in range(1, n)}
Kbar = E.gbar(E.K)
assert (Kbar == np.block([[np.eye(m, dtype=int), np.eye(m, dtype=int)], [np.zeros((m, m), int), np.eye(m, dtype=int)]])).all()
spstore = {}
def sk(Sb):
    k = E.skey(Sb); spstore[k] = Sb; return k
Ik = sk(np.eye(2 * m, dtype=np.int64)); Kk = sk(Kbar)

# ---------------------------------------------------------------- 1. frame game
def frame_states(pr):
    st = []
    for op in pr['ops']:
        if op[0] == 'in': st.append({(op[1], Ik): Fr(1)})
        elif op[0] == 'mix':
            c = {}
            for x, p, r in op[1]:
                Pm = conv(p); assert E.in_P(Pm)
                for (l, g), v in st[r].items():
                    kk = (l, sk(Pm @ spstore[g] % 2)); c[kk] = c.get(kk, 0) + Fr(x) * v
            st.append({k: v for k, v in c.items() if v != 0})
        else:
            z, r = op[1], op[2]
            st.append({(l, sk(Tz[z] @ spstore[g] % 2)): v for (l, g), v in st[r].items()})
    return st
def frame_ok(pr, st):
    for l, rd in enumerate(pr['read']):
        c = {}
        for r, p, x in rd:
            Pm = conv(p); assert E.in_P(Pm)
            for (l2, g), v in st[r].items():
                kk = (l2, sk(Pm @ spstore[g] % 2)); c[kk] = c.get(kk, 0) + Fr(x) * v
        c = {k: v for k, v in c.items() if v != 0}
        if c != {(l, Kk): Fr(1)}: return False
    return True

# ---------------------------------------------------------------- 2. lift to an explicit engine word
Dcache = {}
def Dof(g, mode):
    k = (E.key(g), mode)
    if k not in Dcache: Dcache[k] = E.D(g, mode)
    return Dcache[k]
def Lam_terms(Pm_rows, x, mode):
    """terms (coef, matrix, legit) of x * Lambda(p) = x [phat] D(phat)"""
    Pm = conv(Pm_rows); M, word = E.Ptab[E.skey(Pm)]
    assert (E.gbar(M) == Pm).all()
    out = []
    for (u, j), c in Dof(M, mode).items(): out.append((c.scale(Fr(x)), M @ E.nmat(u, j)))
    return out
def lift(pr, mode):
    steps = []; reg = {}
    for i, op in enumerate(pr['ops']):
        if op[0] == 'in': steps.append(('in', op[1]))
        elif op[0] == 'mix':
            terms = []
            for x, p, r in op[1]:
                terms += [(c, M, reg[r]) for c, M in Lam_terms(p, x, mode)]
            steps.append(('free', terms))
        else:
            z, r = op[1], op[2]
            steps.append(('free', [(c, E.nmat(u, j), reg[r]) for (u, j), c in Dof(E.dirm(z), mode).items()]))
            steps.append(('paid', z, len(steps) - 1))
        reg[i] = len(steps) - 1
    outs = []
    Kinv = E.K.conj().T
    Q = Dof(Kinv, mode) if mode == 'full' else Dof(Kinv, 'f1')
    for l, rd in enumerate(pr['read']):
        terms = []
        for r, p, x in rd: terms += [(c, M, reg[r]) for c, M in Lam_terms(p, x, mode)]
        steps.append(('free', terms)); a = len(steps) - 1
        steps.append(('free', [(c, E.nmat(u, j), a) for (u, j), c in Q.items()])); outs.append(len(steps) - 1)
    return steps, outs

def legit_and_cost(pr, steps):
    paid = [s[1] for s in steps if s[0] == 'paid']
    assert paid == [op[1] for op in pr['ops'] if op[0] == 'move'], "paid letters differ"
    nfree = 0
    for s in steps:
        if s[0] != 'free': continue
        for c, M, src in s[1]:
            Sb = E.gbar(M); assert E.in_P(Sb), "free element outside B"
            nfree += 1
    return len(paid), nfree

# ---------------------------------------------------------------- 3. formal evaluation in C[G], exact
def formal(steps, outs, W):
    regs = []
    for s in steps:
        if s[0] == 'in':
            regs.append([({E.key(E.I): GQ(1)} if l == s[1] else {}) for l in range(W)])
        elif s[0] == 'paid':
            g = E.dirm(s[1]); src = regs[s[2]]
            regs.append([{E.key(g @ E.store[k]): v for k, v in src[l].items()} for l in range(W)])
        else:
            new = [{} for _ in range(W)]
            for c, M, src in s[1]:
                for l in range(W):
                    for k, v in regs[src][l].items():
                        kk = E.key(M @ E.store[k]); new[l][kk] = new[l].get(kk, GQ(0)) + c * v
            regs.append([{k: v for k, v in d.items() if not v.iszero()} for d in new])
    Kk_ = E.key(E.K); ok = True
    for l, o in enumerate(outs):
        for l2 in range(W):
            d = regs[o][l2]
            want = {Kk_: (1, 0)} if l2 == l else {}
            got = {k: (v.r, v.i) for k, v in d.items()}
            if got != want: ok = False
    return ok

# ---------------------------------------------------------------- 4. raw arrays
def apply(M, v, f):
    t = v.reshape((n,) * f)
    for ax in range(f):
        t = np.moveaxis(np.tensordot(M, t, axes=([1], [ax])), 0, ax)
    return t.reshape(-1)
def raw(steps, outs, W, f, seed):
    rn = np.random.default_rng(seed); N = n ** f
    xs = [rn.normal(size=N) + 1j * rn.normal(size=N) for _ in range(W)]
    regs = []
    for s in steps:
        if s[0] == 'in': regs.append(xs[s[1]])
        elif s[0] == 'paid': regs.append(apply(E.dirm(s[1]), regs[s[2]], f))
        else:
            grp = {}
            for c, M, src in s[1]:
                k = (E.key(M), src)
                if k in grp: grp[k][0] += c.c()
                else: grp[k] = [c.c(), M, src]
            out = np.zeros(N, complex)
            for cc, M, src in grp.values():
                if abs(cc) > 1e-15: out += cc * apply(M, regs[src], f)
            regs.append(out)
    err = 0
    for l, o in enumerate(outs):
        t = apply(E.K, xs[l], f); err = max(err, np.abs(regs[o] - t).max() / np.abs(t).max())
    return err

for pr in progs:
    t0 = time.time()
    st = frame_states(pr); fok = frame_ok(pr, st)
    line = f"m={m} {pr['name']:9s} W={pr['W']} M={sum(1 for o in pr['ops'] if o[0]=='move')} frame_identity={fok}"
    for mode in modes:
        steps, outs = lift(pr, mode)
        paid, nfree = legit_and_cost(pr, steps)
        fo = formal(steps, outs, pr['W']) if (m <= 3) else None
        errs = [raw(steps, outs, pr['W'], f, 100 + f) for f in range(1, fmax + 1)]
        line += f" | {mode}: paid={paid} free_terms={nfree} formal_CG={fo} raw_err(f=1..{fmax})=" + \
                ",".join(f"{e:.1e}" for e in errs)
    print(line, f"[{time.time()-t0:.0f}s]", flush=True)
