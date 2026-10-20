#!/usr/bin/env python3
"""orb.py (pk2-orbits, packing round 2, 15 Oct 2026): ORBIT designs in the frame game, with exact screens.

Built on fp-build's fpa.py (frame game, exact end solve, psi screens, constructive charts; copied unchanged into
../frame-game/) and its trunk skeletons (design.py / design_w.py, in ../frame-game/).

An ORBIT injection is a whole sigma-orbit of injections at one site:
    register reg, before its k-th move, add   sum_{i < n} c_i * sigma^i * mu0 * (stored state src)
with sigma in P of order n (n <= NMAX) and c_i independent symbols.  mu0 is the landing chart (fp-build Lemma 3.2)
that puts one walker of src exactly on the frame of a chosen junk walker at a chosen later output of reg.
Kinds of sigma (how the orbit closes):
  stab  : sigma fixes the landing frame Y = w^-1 L0 (all copies land on the junk frame; they differ by tau^i, tau in P)
  cycle : sigma maps the landing frame Y to the landing frame Y' of a SECOND junk walker reachable from the same
          site (wave-7 cycle orbits: the junk frames lie on one sigma-cycle), sigma = chart_to(Y, Y') * stab(Y)
  twof  : sigma in the stabiliser of both L0 and L1 (P cap K^-1 P K, ~ GL(m,2)); the order-3 control is of this kind
  perm  : sigma an address map by a coordinate permutation (commutes with K; o3 type)
  rand  : sigma a random element of P of order <= NMAX
Each design: skeleton + up to D orbit injections (each aimed at an orphan / junk row of the current state).  Then, per
output, three exact verdicts:
  psilin : fp-build's linearised psi screen (every (state, monomial) its own column): FAIL => dead for every scalar
  chlin  : the same relaxation at chart level (end_solve over the monomial components of every state): FAIL => dead
           for every scalar choice, exactly; PASS is necessary only
  exact  : end_solve at drawn rational scalars: YES => a frame-game program (re-verified term by term)
plus heightfilter (p6-prove Theorem C, through a shim for fpa positions): a height system => the word cannot save.
Usage: orb.py MODE m trials depth seed [kinds]     MODE in {two, two2, three}; kinds comma list (default all)
"""
import sys, os, random, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "heights")); sys.path.insert(0, os.path.join(HERE, "..", "frame-game"))
from fpa import *
import design as D
import design_w as DW
import heightfilter as HF

NMAX = 8

# ------------------------------------------------------------------ heightfilter shim (fpa positions as frames)
class FShim:
    def __init__(s, S):
        s.S = S; s.m = S.m; s.L0 = S.pos(S.I); s.LT = S.pos(S.K); s._d = {}
    def d(s, X, Y):
        k = (X, Y)
        if k not in s._d: s._d[k] = s.S.rank_rows(list(X) + list(Y)) - s.S.m
        return s._d[k]
def move_outputs(S, states, W):
    return [{(l, S.pos(g)) for (l, g) in st} for st in states[W:]]
def height_verdict(S, F, states, W):
    outs = move_outputs(S, states, W)
    ok, h = HF.feasible(F, outs, W)
    bb = HF.best_bound(F, outs, W)
    return ok, bb

# ------------------------------------------------------------------ group helpers
def order(S, g, cap=NMAX):
    h = g
    for n in range(1, cap + 1):
        if h == S.I: return n
        h = S.mul(g, h)
    return None
def stab_std(S, d, rng):
    """random element of P fixing U_d = <(e_1|0)..(e_d|0), (0|e_d+1)..(0|e_m)>: addr(A) with A preserving
    span(e_1..e_d), then a phase S_sym with zero top-left d x d block"""
    m = S.m
    while True:
        cols = []
        for k in range(m):
            if k < d: cols.append(rng.randrange(1, 1 << d))           # image of e_k inside span(e_1..e_d)
            else: cols.append(rng.randrange(1, 1 << m))
        if S.rank_rows(cols) == m: break
    A = S.addr(cols)
    Ssym = [0] * m
    for i in range(m):
        for j in range(i, m):
            if i < d and j < d: continue
            if rng.random() < 0.5: Ssym[i] |= 1 << j; Ssym[j] |= 1 << i
    g = S.mul(S.ph(Ssym), A)
    return g
def stab(S, Y, rng):
    """random element of P fixing the frame Y (tuple of basis vectors)"""
    p0, d = normal_P(S, list(Y))
    g = stab_std(S, d, rng)
    h = S.mul(S.inv(p0), S.mul(g, p0))
    assert S.in_P(h) and S._rref([S.apply(h, v) for v in Y]) == S._rref(list(Y))
    return h
def twoframe(S, rng):
    """random element of P cap K^-1 P K: [[A, 0], [A + A^-T, A^-T]]"""
    m = S.m
    while True:
        cols = [rng.randrange(1, 1 << m) for _ in range(m)]
        if S.rank_rows(cols) == m: break
    A = S.addr(cols)                         # [[A,0],[0,A^-T]] (row masks)
    # build [[A,0],[C,A^-T]] with C = A + A^-T: rows m..2m-1 get extra x-bits C
    Arows = [A[i] & ((1 << m) - 1) for i in range(m)]
    AinvT = [(A[m + i] >> m) for i in range(m)]
    rows = list(A)
    for i in range(m):
        rows[m + i] = (A[m + i]) | (Arows[i] ^ AinvT[i])
    g = tuple(rows)
    assert S.is_symp(g) and S.in_P(g) and S.in_P(S.mul(S.K, S.mul(g, S.inv(S.K))))
    return g
def perm_elt(S, rng):
    m = S.m; p = list(range(m)); rng.shuffle(p)
    return S.addr([1 << p[k] for k in range(m)])

# ------------------------------------------------------------------ split a symbolic state into monomial components
def mono_split(states):
    out = []
    for st in states:
        bym = {}
        for k, v in st.items():
            pv = v if isinstance(v, Poly) else Poly.const(v)
            for mono, co in pv.t.items():
                d = bym.setdefault(mono, {}); d[k] = d.get(k, 0) + co
        for mono, d in bym.items():
            d = {k: v for k, v in d.items() if v != 0}
            if d: out.append(d)
    return out
def subs(states, vals):
    out = []
    for st in states:
        d = {}
        for k, v in st.items():
            if isinstance(v, Poly):
                x = Fr(0)
                for mono, co in v.t.items():
                    t = co
                    for (i, e) in mono: t *= vals[i] ** e
                    x += t
            else: x = v
            if x != 0: d[k] = x
        out.append(d)
    return out

# ------------------------------------------------------------------ orbit designs
class Orbit:
    def __init__(s, S, P, W, rng, kinds):
        s.S, s.P, s.W, s.rng, s.kinds = S, P, W, rng, kinds
        s.nvar = 0
        s.log = []
    def sites(s, l, X):
        """(reg, k, t, src, g, w) with walker g (input l) of stored state src landable on frame X at output t of reg"""
        S, P = s.S, s.P
        out = []
        for reg in P.order:
            n = len(P.regs[reg])
            for k in range(n):
                w = S.I
                for t in range(k, n):
                    w = S.mul(S.transv(P.regs[reg][t]), w)
                    lw = S.level(w)
                    for src in range(P.avail[(reg, k)]):
                        for (l2, g), v in P.states[src].items():
                            if l2 != l: continue
                            if S.rank_rows(list(S.pos(g)) + list(X)) - S.m != lw: continue
                            out.append((reg, k, t, src, g, w))
        return out
    def add_orbit(s, row, others):
        """aim an orbit at junk row (l, element); others = other junk rows (for cycle orbits)"""
        S, rng = s.S, s.rng
        l, j = row; X = S.pos(j)
        cands = s.sites(l, X)
        rng.shuffle(cands)
        for (reg, k, t, src, g, w) in cands[:60]:
            mu0 = find_chart_c(S, w, g, X)
            if mu0 is None: continue
            Y = S._rref([S.apply(S.inv(w), 1 << (S.m + i)) for i in range(S.m)])   # w^-1 L0
            if s.kinds == ["cyclefirst"]:
                # pure cycle orbits need two junk walkers of one input: one stabiliser orbit first, then cycles only
                kind = "stab" if not s.log else "cycle"
            else:
                kind = rng.choice(s.kinds)
            sig = None
            for _ in range(40):
                if kind == "stab":
                    # sigma fixes the frame mu0 g X... in the chart's coordinates: the copies sigma^i mu0 g land on X
                    # iff sigma fixes Y = w^-1 L0 (then pos(w sigma^i mu0 g) = g^-1 mu0^-1 sigma^-i Y = X)
                    c = stab(S, Y, rng)
                elif kind == "cycle":
                    # second junk frame X' of the same input: the frame Y' = mu0 g X' (so copy 1 lands on X')
                    oth = [r for r in others if r[0] == l and S.pos(r[1]) != X]
                    if not oth: c = None; break
                    j2 = rng.choice(oth)[1]; X2 = S.pos(j2)
                    Y2 = S._rref([S.apply(S.mul(mu0, g), v) for v in X2])
                    # want sigma^-1 Y = Y2  (copy 1: pos = g^-1 mu0^-1 sigma^-1 Y = X2)
                    if S.rank_rows(list(Y2) + list(S.pos(S.I))) != S.rank_rows(list(Y) + list(S.pos(S.I))):
                        c = None; break
                    base = chart_to(S, Y2, Y)        # base(Y2) = Y
                    if base is None: c = None; break
                    c = S.mul(base, stab(S, Y2, rng))   # c(Y2) = Y, i.e. c^-1 Y = Y2
                elif kind == "twof":
                    c = twoframe(S, rng)
                elif kind == "perm":
                    c = perm_elt(S, rng)
                else:
                    c = S.rand_P(rng)
                if c is None: break
                n = order(S, c)
                if n is not None and n >= 2: sig = (c, n); break
            if sig is None: continue
            c, n = sig
            h = S.I
            for i in range(n):
                s.P.inj.append((reg, k, src, S.mul(h, mu0), Poly.var(s.nvar))); s.nvar += 1
                h = S.mul(c, h)
            s.log.append((kind, reg, k, t, src, n))
            return True
        return False

def junk_rows(S, states, out_idx, W):
    """junk walkers of the delivering states: every walker except (l, K-frame walker) of state out_idx[l]"""
    rows = []
    for l in range(W):
        st = states[out_idx[l]]
        for (l2, g), v in st.items():
            if l2 == l and S.pos(g) == S.pos(S.K): continue
            rows.append((l2, g))
    return rows

CHAR = "char" in sys.argv
def draw_vals(rng, nvar, olog):
    """generic rationals, or (with 'char') each orbit's coefficients = overall rational x a rational character of
    <sigma> (trivial, or the sign character for even orbit length): the closed-pattern coefficient choices"""
    if not CHAR or not olog:
        return [Fr(rng.choice([1, -1, 2, -2, 3, -3, 5])) / rng.choice([1, 2, 3, 7]) for _ in range(max(nvar, 1))]
    vals = []
    for o in olog:
        n = o[5]; c = Fr(rng.choice([1, -1, 2, -2, 3])) / rng.choice([1, 2, 3])
        sign = (n % 2 == 0) and rng.random() < 0.5
        vals += [c * (-1 if (sign and i % 2) else 1) for i in range(n)]
    return vals + [Fr(1)] * max(0, max(nvar, 1) - len(vals))

def evaluate(S, F, P, W, nvar, rng, do_exact=True, olog=None):
    st = P.run()
    res = {}
    res["psilin"] = "".join("1" if psi_lin_screen(S, st, l) else "0" for l in range(W))
    split = mono_split(st)
    ch = []
    for l in range(W):
        r = end_solve(S, split, {(l, S.K): Fr(1)}, cap_cols=6000, cap_rows=20000)
        ch.append(r[0][0])
    res["chlin"] = "".join(ch)
    hv, bb = height_verdict(S, F, st, W)
    res["height"] = "H" if hv else "-"
    res["hbound"] = bb
    if do_exact:
        vals = draw_vals(rng, nvar, olog)
        stn = subs(st, vals)
        ex = []
        for l in range(W):
            r = end_solve(S, stn, {(l, S.K): Fr(1)}, cap_cols=6000, cap_rows=20000)
            if r[0] == "yes":
                ok, diff = verify_combo(S, stn, r[1]["solution"], {(l, S.K): Fr(1)})
                ex.append("Y" if ok else "?")
            else:
                ex.append(r[0][0])
        res["exact"] = "".join(ex)
        res["vals"] = vals
        # psi at the DRAWN scalars (exact; respects rigidity, unlike psilin): necessary for these scalars
        res["psinum"] = "".join("1" if (ex[l] == "Y" or psi_screen(S, stn, l)) else "0" for l in range(W))
        # at least W states must have walkers at level m (dimension count, PROOF.md Remark 7); else skip the layers
        ntop = sum(1 for x in stn if any(S.level(g) == S.m for (_, g) in x))
        res["ntop"] = ntop
        if ntop < W:
            res["top"] = "T%d" % ntop
        else:
            res["top"] = ",".join(layer_screen(S, stn, l, cap_cols=8000, cap_rows=40000) if ex[l] == "c" else "-"
                                  for l in range(W))
    return res, st

def top_screen(S, states, l, cap_cols=40000, cap_rows=120000):
    """EXACT necessary condition (any scalars as given): K e_l lies in the C[P]-span of the level-m parts of the
    stored states (left P-multiplication keeps levels, so the readout's level-m part must already produce K e_l)."""
    m = S.m
    top = [{k: v for k, v in st.items() if S.level(k[1]) == m} for st in states]
    top = [t for t in top if t]
    if not top: return "n"
    r = end_solve(S, top, {(l, S.K): Fr(1)}, cap_cols=cap_cols, cap_rows=cap_rows)
    return r[0][0] if r[0] != "yes" else "y"

def layer_screen(S, states, l, cap_cols=20000, cap_rows=80000):
    """Layered EXACT necessary conditions (wave-7 style, here in the frame game): for d = m, m-1, ..., 0 the readout
    restricted to walkers of level >= d must already give K e_l (left P-multiplication keeps levels).  Returns a
    string: 'L<d>n' = first layer d that FAILS (so the design is exact NO at these scalars), 'L<d>c' = layer d over
    the cap (layers above passed), 'L0y' = every layer passed (then the full solve decides)."""
    for d in range(S.m, -1, -1):
        sub = [{k: v for k, v in st.items() if S.level(k[1]) >= d} for st in states]
        sub = [t for t in sub if t]
        if not sub: return "L%dn" % d
        r = end_solve(S, sub, {(l, S.K): Fr(1)}, cap_cols=cap_cols, cap_rows=cap_rows)
        if r[0] == "no": return "L%dn" % d
        if r[0] == "cap": return "L%dc" % d
    return "L0y"

def make_prog(S, m, W, rng, mode):
    if mode == "three":
        P = DW.ProgW(S, m, W, rng)
    else:
        P = D.Prog(S, m, rng, dup=(rng.randrange(0, m) if rng.random() < 0.5 else 0), tlen=(2 if mode == "two2" else 1))
    return P

def deliver_idx(P, W):
    regs = ["A", "B"] if W == 2 else ["A%d" % l for l in range(W)]
    out = []
    for reg in regs:
        n = len(P.regs[reg])
        out.append(P.out[(reg, n - 1)] if n else P.out[("Rt", len(P.regs["Rt"]) - 1)])
    return out

def main():
    mode, m, trials, depth, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    kinds = sys.argv[6].split(",") if len(sys.argv) > 6 and sys.argv[6] not in ("all", "char") else ["stab", "cycle", "twof", "perm", "rand"]
    W = 3 if mode == "three" else 2
    S = Sp(m); F = FShim(S); rng = random.Random(seed)
    tl = 2 if mode == "two2" else 1
    M = W * (m - 1) + 1 if mode == "three" else 2 * m - tl
    print(f"orb.py mode={mode} m={m} W={W} target {W} in {M}; trials={trials} depth={depth} seed={seed} kinds={kinds}",
          flush=True)
    tally = {}; t0 = time.time(); nfound = 0
    for tr in range(trials):
        P = make_prog(S, m, W, rng, mode)
        if not P.ok: tally["skel-fail"] = tally.get("skel-fail", 0) + 1; continue
        O = Orbit(S, P, W, rng, kinds)
        hist = []
        for step in range(depth + 1):
            res, st = evaluate(S, F, P, W, O.nvar, rng, do_exact=True, olog=O.log)
            hist.append((res["psilin"], res["chlin"], res["exact"], res["height"], res["hbound"], res["top"],
                         res["psinum"]))
            if all(x == "Y" for x in res["exact"]):
                nfound += 1
                print("FRAME-GAME PROGRAM FOUND", mode, m, "trial", tr, "log", O.log, "vals", res["vals"], flush=True)
                import pickle
                with open(os.path.join(HERE, f"found_{mode}_m{m}_s{seed}_t{tr}.pkl"), "wb") as fh:
                    pickle.dump(dict(m=m, mode=mode, regs=P.regs, order=P.order, inj=P.inj, vals=res["vals"],
                                     a=getattr(P, "a", (P.a1, P.a2)), p2=getattr(P, "p2", None),
                                     p=getattr(P, "p", None)), fh)
                break
            if step == depth: break
            st_num = subs(st, res["vals"])
            dl = deliver_idx(P, W)
            rows = junk_rows(S, st_num, dl, W)
            # prefer orphan rows of the exact closure (as fp-build), else the delivering states' junk
            orph = []
            for l in range(W):
                r = end_solve(S, st_num, {(l, S.K): Fr(1)}, cap_cols=6000, cap_rows=20000)
                if r[0] == "no" and "singles" in r[1] and r[1]["singles"]:
                    dmin = r[1]["singles"][0][0]
                    orph += [row for (dd, row) in r[1]["singles"] if dd <= dmin + 1]
            target_rows = orph if orph else rows
            if not target_rows: break
            rng.shuffle(target_rows)
            done = False
            for row in target_rows[:4]:
                if O.add_orbit(row, rows + orph): done = True; break
            if not done: break
        last = hist[-1]
        key = (mode, "a=" + str(getattr(P, "a", (P.a1, P.a2))), "orbits=%d" % len(O.log),
               "psilin=" + last[0], "chlin=" + last[1], "exact=" + last[2], "height=" + last[3])
        tally[key] = tally.get(key, 0) + 1
        kk = ("KIND", tuple(sorted(set(x[0] for x in O.log))), "psilin=" + last[0], "chlin=" + last[1])
        tally[kk] = tally.get(kk, 0) + 1
        print(f"trial {tr}: a={getattr(P, 'a', (P.a1, P.a2))} dup={getattr(P, 'dup', 0)} orbits={O.log} "
              f"hist={hist}", flush=True)
    print(f"done: {trials} trials in {time.time() - t0:.0f}s; found {nfound}", flush=True)
    for k in sorted(tally, key=str): print("  ", k, tally[k], flush=True)

if __name__ == "__main__":
    main()
