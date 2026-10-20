"""lemma5_check.py -- w3-asym: exact check of Lemma 5 on random instances, OUTPUT 0 ALONE, merge coefficient forced
nonzero (Rabinowitsch: a0*z = 1).  (A) cores with lambda(J0) = 3 plus 1-3 random injections of every allowed kind
(x_j -> u; x_j, s_j -> v_e; x_j, s_j, u -> t_e): Lemma 5 predicts INCONSISTENT.  (B) controls: cores with lambda(J0) = 2
plus the landing injection x_1 -> u read through v_0 (frame of output 0's junk): output 0 alone may be CONSISTENT.
Writes Singular scripts.  Usage: lemma5_check.py nA nB seed cap outdir"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_decide import *

def core_lam(rng, want, W=2, tries=3000):
    for _ in range(tries):
        bs = [prog_letter(rng, [G.ID]) for l in range(W)]
        hs = [G.intern(G.mul(G.D1, b)) for b in bs]
        as_ = [prog_letter(rng, [h]) for h in hs]
        ps = [G.mul(G.D1, a, G.MATS[h]) for a, h in zip(as_, hs)]
        if lam(G.mul(G.K, np.conj(ps[0]).T, ps[1])) == want:
            break
    else:
        return None
    pids = [G.intern(p) for p in ps]
    rec = [[(l, bs[l])] for l in range(W)]
    rec.append([(W + l, as_[l]) for l in range(W)])
    u = 2 * W
    for l in range(W):
        path = midpaths(rng, pids[l], nsamp=300)
        if not path: return None
        rec.append([(u, path[0][0])]); rec.append([(u + 1 + 2 * l, path[0][1])])
    return rec, ps

def script_out0(P, cap):
    polys = monomial_split(P)
    names = [f"a{k}" for k in range(P.np)]
    C = Closure(P, 0, cap)
    if not C.ok: return None, None
    rows = {}
    for (r, c, i, j) in C.E: rows.setdefault(r, []).append((c, i, j))
    bn = [f"b{c}" for c in range(C.nc)]
    eqs = []
    for r in range(C.nr):
        parts = [f"({poly_str(polys[i], j, names)})*{bn[c]}" for (c, i, j) in rows.get(r, [])]
        e = "+".join(parts) if parts else "0"
        if r == 0: e += "-1"
        eqs.append(e)
    eqs.append("a0*z-1")
    s = f"ring R = 0, ({','.join(bn + names + ['z'])}), dp;\nideal I = " + ",\n  ".join(eqs) + ";\n"
    s += "option(redSB);\nideal J = std(I);\nif (size(J) == 1 and J[1] == 1) { print(\"INCONSISTENT\"); } else { print(\"CONSISTENT\"); }\nquit;\n"
    return s, (C.nr, C.nc)

if __name__ == '__main__':
    nA, nB, seed, cap, outdir = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed); W = 2
    k = 0
    while k < nA:
        r = core_lam(rng, 3)
        if r is None: continue
        rec, ps = r
        r2 = [list(m) for m in rec]
        u = 2 * W; vmoves = [W + 1, W + 3]; tmoves = [W + 2, W + 4]
        for _ in range(rng.randint(1, 3)):
            kind = rng.choice(['xu', 'xv', 'sv', 'xt', 'st', 'ut'])
            j = rng.randrange(W); e = rng.randrange(W)
            if kind == 'xu': r2[W].append((j, rand_engine(rng)))
            elif kind == 'xv': r2[vmoves[e]].append((j, rand_engine(rng)))
            elif kind == 'sv': r2[vmoves[e]].append((W + j, rand_engine(rng)))
            elif kind == 'xt': r2[tmoves[e]].append((j, rand_engine(rng)))
            elif kind == 'st': r2[tmoves[e]].append((W + j, rand_engine(rng)))
            else: r2[tmoves[e]].append((u, rand_engine(rng)))
        P = build_rec(W, r2)
        scr, sz = script_out0(P, cap)
        if scr is None: continue
        open(os.path.join(outdir, f"A_{seed}_{k}.sing"), 'w').write(scr); print(f"A {k}: params {P.np} closure {sz}", flush=True); k += 1
    k = 0
    while k < nB:
        r = core_lam(rng, 2)
        if r is None: continue
        rec, ps = r
        J = G.mul(G.K, np.conj(ps[0]).T, ps[1]); fk = G.frame_key(J)
        ks = landing_letters(rng, fk, G.mul(G.D1, rec[W + 1][0][1]), np.eye(16, dtype=complex), 3000)
        if not ks: continue
        r2 = [list(m) for m in rec]; r2[W].append((1, ks[0]))
        P = build_rec(W, r2)
        scr, sz = script_out0(P, cap)
        if scr is None: continue
        open(os.path.join(outdir, f"B_{seed}_{k}.sing"), 'w').write(scr); print(f"B {k}: params {P.np} closure {sz}", flush=True); k += 1
