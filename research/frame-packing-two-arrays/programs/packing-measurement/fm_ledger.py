#!/usr/bin/env python3
"""fp-measure: the raw-gauge ledger (sl-combinatorics PROOF.md Definition 1.4, Theorem 2.1) of one call of a local
word in gcert/1 format (the hub-line words U8X, X4X, a T74-line word, our E8 unit r783), in the frame model.

Frame model (sl-combinatorics Lemma 2.5; derivation in NOTE.md section 2):
- a register at subspace frame U sits at the Lagrangian F_U = {(c | a + c) : c in U, a in U^perp};
- an input whose source frame (twist) is S has, in a register at frame U, the raw label Phi_U Phi_S^-1, whose level
  is d(F_U, F_S) = dim U + dim S - 2 dim(U cap S);
- two inputs share a raw label in a register exactly when they have the same source frame S (same twist);
- a nested climb U -> U' of k paid moves passes through nested frames only (geodesic betweenness), so each label
  moves +-1 per move (Z = 0), and label S descends exactly dim(U' cap S) - dim(U cap S) times; a descent U -> U'
  (U' inside U) gives dim(U cap S) - dim(U' cap S) ascents.
Source frames in the window: x_t on its port line <p_t> (distinct per t); y_t at 0 (ONE class: identity twist);
slot dirt (the helper's own data, when hosted) at 0, either one class per slot (dirt='*-each') or one class for all
(dirt='*-one'), on the exterior-gauged slots only (dirt='ext-*': the hosted helpers, whose dirt E7 keeps away from y)
or on every slot (dirt='all-*', a reference: then dirt reaches y), or no dirt at all (dirt='none': the gcert model,
slots start empty).

Centre model (an assumption, flagged in NOTE.md): each retained total k at frame f is copied (free) into a fresh
centre register, which DESCENDS f -> 0 in dim f paid moves (this is the cst = sum dim f that the local price charges),
and is then scattered into the targets at frame 0 with the table's coefficients; the centre is left as scratch.

Exact rational contents (Fractions).  The segment identity M = Psi_end - Psi_0 + Bmix + 2D + Z - P + I0 is asserted.
usage: python3 -I fm_ledger.py WORD.gcert.json.gz NAME [dirt]   (dirt in none, ext-each, ext-one, all-each, all-one)"""
import sys, json, gzip
from fractions import Fraction as Fr
from collections import Counter, defaultdict


def load(p):
    return json.load(gzip.open(p, "rt") if p.endswith(".gz") else open(p))


def span_has(basis, u):
    """u in span(basis); basis is reduced echelon (distinct pivots), as in gcert frames"""
    for b in basis:
        if u >> (b.bit_length() - 1) & 1:
            u ^= b
    return u == 0


def run(D, dirt="none", centre="descend", init="F1", pos=False):
    h, v, R = D["h"], D["v"], D["R"]
    fr = D["frames"]
    dim = [len(f) for f in fr]
    ports = D["ports"]
    n = 2 * v + R
    nreg = n + len(D["ret"])                       # extra registers: the centre copies
    cur = list(D["start"]) + [None] * len(D["ret"])
    final = list(D["final"])
    # contents: register -> {input key: Fraction}; input keys ('x',t), ('y',t), ('d',r)
    reg = [dict() for _ in range(nreg)]
    for t in range(v):
        if init == "F1":      # fresh forward call: x_t holds x_t, y_t holds y_t
            reg[t] = {("x", t): Fr(1)}
            reg[v + t] = {("y", t): Fr(1)}
        elif init == "F2":    # forward call after a backward one: x_t holds x_t - y_t, y_t holds y_t
            reg[t] = {("x", t): Fr(1), ("y", t): Fr(-1)}
            reg[v + t] = {("y", t): Fr(1)}
        elif init == "F3":    # third shear of a swap: x_t holds -y_t, y_t holds y_t + x_t
            reg[t] = {("y", t): Fr(-1)}
            reg[v + t] = {("y", t): Fr(1), ("x", t): Fr(1)}
        else:
            raise ValueError(init)
    x0 = [dict(reg[t]) for t in range(v)]
    y0 = [dict(reg[v + t]) for t in range(v)]
    extset = set(r for r, _ in D.get("ext", []))
    if dirt != "none":
        for r in range(2 * v, n):
            if dirt.startswith("all") or dirt.startswith("gates") or r in extset:
                reg[r] = {("d", r): Fr(1)}
    # class of an input key, and its source frame as (dim S, line vector or None)
    def cls(key):
        k, i = key
        if k == "x":
            return ("x", i)
        if k == "y":
            return ("Y",)
        return ("d", i) if (dirt.endswith("each") or dirt == "gates") else ("D",)   # gates-one: one dirt class
    memb = {}
    # positions (vertices of Gamma_h) of the raw labels: Phi_S F_U, F_U = <(c|c), (0|a)> for c in U, a in U^perp;
    # vectors (a|b) as ints with the x part in the high h bits; Phi_<p> acts as T_p: (a|b) -> (a + (p.b) p | b)
    perpc, Fc, posc = {}, {}, {}
    def perp(f):
        if f not in perpc:
            perpc[f] = [a for a in range(1, 1 << h) if all(bin(a & c).count("1") % 2 == 0 for c in fr[f])]
        return perpc[f]
    def echelon(vs):
        piv = {}
        for x in vs:
            for pb in sorted(piv, reverse=True):
                if x >> pb & 1:
                    x ^= piv[pb]
            if x:
                pb = x.bit_length() - 1
                for q in list(piv):
                    if piv[q] >> pb & 1:
                        piv[q] ^= x
                piv[pb] = x
        return tuple(sorted(piv.values()))
    def FU(f):
        if f not in Fc:
            vs = [(c << h) | c for c in fr[f]]
            pb = echelon(perp(f))
            vs += list(pb)
            Fc[f] = echelon(vs)
        return Fc[f]
    def position(c, f):
        kk = (c, f)
        if kk not in posc:
            if c[0] != "x":
                posc[kk] = FU(f)
            else:
                p = ports[c[1]]
                out = []
                for w in FU(f):
                    a, b = w >> h, w & ((1 << h) - 1)
                    if bin(p & b).count("1") & 1:
                        a ^= p
                    out.append((a << h) | b)
                posc[kk] = echelon(out)
        return posc[kk]
    poshist = Counter()
    Pp = [0]
    def inU(p, f):
        kk = (p, f)
        if kk not in memb:
            memb[kk] = span_has(fr[f], p)
        return memb[kk]
    def level(c, f):
        """d(F_U, F_S), U = frame f"""
        if c[0] == "x":
            return dim[f] + 1 - 2 * (1 if inU(ports[c[1]], f) else 0)
        return dim[f]
    def classes(r):
        return set(cls(k) for k in reg[r])
    def psi(r):
        f = cur[r]
        return sum(level(c, f) for c in classes(r)) if reg[r] else 0
    L = dict(M=0, A=0, D=0, Z=0, P=0, I0=0, Bmix=Fr(0))
    by = defaultdict(Counter)                      # per register kind: M, P, D, I0
    hist = Counter()                               # number of labels on a move -> moves
    maxlab = 0
    movelog = []                                   # (register, classes on the move, number of moves)
    def kind(r):
        return "x" if r < v else "y" if r < 2 * v else "s" if r < n else "c"
    def move(r, f_new, down=False):
        nonlocal maxlab
        f_old = cur[r]
        if f_old == f_new:
            return
        k = abs(dim[f_new] - dim[f_old])
        if not down:
            assert dim[f_new] > dim[f_old] and all(span_has(fr[f_new], u) for u in fr[f_old]), "not a nested climb"
        cl = classes(r)
        L["M"] += k
        by[kind(r)]["M"] += k
        if not cl:
            L["I0"] += k
            by[kind(r)]["I0"] += k
            hist[0] += k
        else:
            L["P"] += k * (len(cl) - 1)
            by[kind(r)]["P"] += k * (len(cl) - 1)
            hist[len(cl)] += k
            if pos:
                nonlocal_pos = len(set(position(c, f_old) for c in cl))
                poshist[nonlocal_pos] += k
                Pp[0] += k * (nonlocal_pos - 1)
            maxlab = max(maxlab, len(cl))
            movelog.append((r, frozenset(cl), k))
            for c in cl:
                if c[0] == "x":
                    a, b = (1 if inU(ports[c[1]], f_old) else 0), (1 if inU(ports[c[1]], f_new) else 0)
                    dd = (b - a) if not down else (k - (a - b))   # descents of this label over the k moves
                else:
                    dd = 0 if not down else k
                assert 0 <= dd <= k
                L["D"] += dd
                L["A"] += k - dd
                by[kind(r)]["D"] += dd
        cur[r] = f_new
    def mix(targets, fn):
        before = sum(psi(r) for r in targets)
        fn()
        for r in targets:
            reg[r] = {kk: c for kk, c in reg[r].items() if c != 0}
        L["Bmix"] -= sum(psi(r) for r in targets) - before
    def axpy(t, src, k):
        d = reg[t]
        for kk, u in src.items():
            d[kk] = d.get(kk, Fr(0)) + k * u
    psi0 = sum(psi(r) for r in range(nreg) if cur[r] is not None)
    if dirt.startswith("gates"):
        # the engine's start dirt gate g1: y -= K s at the start frames (K read off a dirty forward replay)
        import os as _os
        sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
        from bw_ledger import forward_matrices
        Kd, _, _ = forward_matrices(D)
        tg = [t for t in range(v) if Kd[t]]
        for t in tg:
            assert cur[v + t] == 0 and all(cur[s_] == 0 for s_ in Kd[t])
        mix([v + t for t in tg], lambda: [axpy(v + t, dict(reg[s_]), -c_) for t in tg for s_, c_ in Kd[t].items()])
    def regs(g):
        return [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == "out" else [])
    def gate(g):
        f = g[1]
        for r in regs(g):
            move(r, f)
        if g[0] == "out":
            src = dict(reg[g[2]])
            tg = [t for t, _, _ in g[3]]
            mix(tg, lambda: [axpy(t, src, Fr(a, b)) for t, a, b in g[3]])
        else:
            mix([g[2]], lambda: [axpy(g[2], dict(reg[s]), Fr(a, b)) for s, a, b in g[3]])
    for g in D["A"]:
        gate(g)
    # centre: copy each retained total into a fresh register at its frame, descend to 0, scatter into y at 0
    cregs = []
    for k, (_, s, f) in enumerate(sorted(D["ret"])):
        assert cur[s] == f
        c = n + k
        cur[c] = f
        mix([c], lambda: axpy(c, dict(reg[s]), Fr(1)))
        move(c, 0, down=True)
        cregs.append(c)
    if "table" in D["scat"]:
        sc = D["scat"]["table"]
    else:                                          # star rule: one total per coordinate (gxcore.scat_row)
        ins, outs = D["scat"]["inside"], D["scat"]["outside"]
        sc = [[(k,) + tuple(ins if ports[t] >> k & 1 else outs) for k in range(h)] for t in range(v)]
    for t in range(v):
        assert cur[v + t] == 0
        mix([v + t], lambda: [axpy(v + t, dict(reg[cregs[k]]), Fr(a, b)) for k, a, b in sc[t]])
    for g in D["B"]:
        gate(g)
    for r in range(n):
        move(r, final[r])
    if dirt.startswith("gates"):
        # the end dirt gate g8: slots := rows of the inverse, which restores every slot to its own dirt
        def g8():
            for r in range(2 * v, n):
                reg[r] = {("d", r): Fr(1)}
        mix(list(range(2 * v, n)), g8)
    psi_end = sum(psi(r) for r in range(nreg))
    # exactness of the local word: x restored, y_t = x_t + y_t
    okx = all(reg[t] == x0[t] for t in range(v))
    def plus(a, b):
        o = dict(a)
        for kk, u in b.items():
            o[kk] = o.get(kk, 0) + u
        return {kk: u for kk, u in o.items() if u != 0}
    oky = all(reg[v + t] == plus(y0[t], x0[t]) for t in range(v))
    # delivered levels: x_t in x_t, x_t and y_t in y_t
    deliv = sum(level(c, final[t]) for t in range(v) for c in set(cls(k) for k in reg[t])) + \
        sum(level(c, final[v + t]) for t in range(v) for c in set(cls(k) for k in reg[v + t]))
    scratch = psi_end - sum(psi(r) for r in range(2 * v))      # what slots and centres carry at the end
    dirt_end = Counter()
    for r in range(2 * v, n):
        for kk in reg[r]:
            dirt_end["slot terms"] += 1
    # kept vs junk: a label on a move is "kept" if the register still holds that class at the end of the call
    fin_cls = [classes(r) for r in range(nreg)]
    P_keep = sum(k * max(0, len(cl & fin_cls[r]) - 1) for r, cl, k in movelog)
    lab_moves = sum(k * len(cl) for r, cl, k in movelog)
    junk_moves = sum(k * len(cl - fin_cls[r]) for r, cl, k in movelog)
    rhs = psi_end - psi0 + L["Bmix"] + 2 * L["D"] + L["Z"] - L["P"] + L["I0"]
    assert rhs == L["M"], ("segment identity fails", rhs, L["M"])
    out = dict(L)
    out["Bmix"] = int(L["Bmix"]) if L["Bmix"].denominator == 1 else str(L["Bmix"])
    out.update(psi0=psi0, psi_end=psi_end, deliv=deliv, scratch_end=scratch, x_restored=okx, y_gets_x=oky,
               identity_ok=True, maxlab=maxlab, hist=dict(sorted(hist.items())),
               by={k: dict(c) for k, c in by.items()}, slot_terms_end=dirt_end["slot terms"],
               # call-level B: level destroyed at mixes + level carried away by scratch and junk beyond the delivery
               B=int(L["Bmix"] + psi_end - deliv - (psi0 if False else 0)),
               net=int(psi_end - psi0 - L["M"]), P_pos=Pp[0] if pos else None,
               pos_hist=dict(sorted(poshist.items())) if pos else None, P_keep=P_keep, label_moves=lab_moves, junk_label_moves=junk_moves)
    out["net_identity"] = "P - Bmix - 2D - Z - I0 = Psi_end - Psi_0 - M"
    return out


def main():
    p, name = sys.argv[1], sys.argv[2]
    modes = sys.argv[3:] or ["none", "ext-each", "ext-one", "all-each"]
    inits = ["F1"]
    if modes and modes[0].startswith("init="):
        inits = modes[0][5:].split(","); modes = modes[1:] or ["none"]
    D = load(p)
    res = {}
    for ini in inits:
      for dm in modes:
        r = run(D, dirt=dm, init=ini)
        dm = dm if ini == "F1" else "%s/%s" % (ini, dm)
        res[dm] = r
        M = r["M"]
        print("%s dirt=%s: M=%d (local N %d, diff %d) P=%d Bmix=%s D=%d Z=%d I0=%d A=%d | Psi0=%d Psi_end=%d deliv=%d"
              " scratch_end=%d | P/M=%.3f D/M=%.3f maxlab=%d | x restored %s, y_t = x_t + y_t %s"
              % (name, dm, M, D["N"], M - D["N"], r["P"], r["Bmix"], r["D"], r["Z"], r["I0"], r["A"], r["psi0"],
                 r["psi_end"], r["deliv"], r["scratch_end"], r["P"] / M, r["D"] / M, r["maxlab"], r["x_restored"],
                 r["y_gets_x"]))
        for k in "xysc":
            b = r["by"].get(k, {})
            if b:
                print("    %s: M=%d P=%d D=%d I0=%d  P/M=%.2f" % (k, b.get("M", 0), b.get("P", 0), b.get("D", 0),
                                                              b.get("I0", 0), b.get("P", 0) / max(1, b.get("M", 0))))
        print("    labels per move (moves weighted):", r["hist"])
        print("    kept packing P_keep=%d (labels the register still holds at the end); junk share of label-moves %.3f"
              % (r["P_keep"], r["junk_label_moves"] / r["label_moves"]))
    json.dump(dict(word=p, name=name, h=D["h"], v=D["v"], R=D["R"], N=D["N"], cst=D["cst"], modes=res),
              open("out/fm_%s.json" % name, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
