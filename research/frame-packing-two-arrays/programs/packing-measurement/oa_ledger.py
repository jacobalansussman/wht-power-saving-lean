#!/usr/bin/env python3
"""fp-measure: the raw-gauge ledger (P, B, D, Z, I0) of OpenAI's complex exchange network (explicit-DFT preprint,
Section 2.1, rows (2.3) and Lemma 2.1) at small ground sets h, in the frame model.

Frames: an-openai's labels_check.py construction (explicit subspaces of F_2^(h^3), one invocation per stage, the
representative target S = the invocation's own triple at the varying position, every T enumerated).  Our own code;
Sp and perp_within are imported from labels_check.py beside this file (our own program, copied unchanged).

Contents (derived in NOTE.md section 4 from rows (2.3) and Lemma 2.1; stage 1 and 3 forward X -> Y, stage 2 the
inverse invocation Y -> X, rows run 7..0 with negated updates):
- raw classes: every X input d has its own twist (source line <u_d>); Y inputs, centre and side wires all have the
  identity twist (source 0), so all of them form ONE class 'Y' in any register;
- the labels a wire holds on each edge between two touches are listed in EDGES below.
Raw level of class S in a register at frame L: dim L + dim S - 2 dim(L cap S) (Lemma 2.5 of sl-combinatorics).
Checks: M = the brute-force move sum of labels_check (W m - 2N + 2L); every wire ends with ONE label at level m
(so B = Bmix); and the unit identity  Wm - M = P - Bmix - 2D - Z - I0, with Bmix computed independently from the
level changes at the gates.
usage: python3 -I oa_ledger.py h [h ...]"""
import sys, os, itertools, importlib.util, json
from math import comb

spec = importlib.util.spec_from_file_location(
    "lc", os.path.join(os.path.dirname(os.path.abspath(__file__)), "labels_check.py"))
lc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lc)
Sp, perp_within, dot = lc.Sp, lc.perp_within, lc.dot


def main(h):
    T = list(itertools.combinations(range(h), 3))
    ind = {S: sum(1 << i for i in S) for S in T}
    v = len(T)
    nbr = lambda S, U: S != U and len(set(S) & set(U)) % 2 == 0
    m = h ** 3
    E1 = [1 << i for i in range(h)]
    a_fix = [T[0], T[5], T[11]]
    tvec = [ind[a] for a in a_fix]

    def flat(fs):
        """tensor of factor vectors (list of h-bit ints) -> int over h^len bits, first factor most significant"""
        out = 1
        nb = 0
        res = 1
        # build by repeated Kronecker
        res, nb = 1, 0
        for f in fs:
            new = 0
            for i in range(nb if nb else 1):
                pass
            # kron(res (nb bits, or scalar), f (h bits))
            if nb == 0:
                res, nb = f, h
                continue
            new = 0
            for i in range(nb):
                if (res >> i) & 1:
                    new |= f << (i * h)
            res, nb = new, nb * h
        return res

    def tens_space(parts):
        """parts: list (one per position 0..2) of factor-space bases (lists of h-bit ints), or ('flatpre', basis, n)
        for a subspace of the first n positions given in flat coordinates.  Returns Sp of the tensor product."""
        # expand: combine consecutive factor lists into flat bases
        bases = []          # list of (basis of flat vectors, number of positions)
        for p in parts:
            if isinstance(p, tuple):
                bases.append((p[1], p[2]))
            else:
                bases.append((p, 1))
        cur, npos = [1], 0
        for bs, k in bases:
            nb_new = h ** k
            new = []
            for a in cur:
                for b in bs:
                    if npos == 0:
                        new.append(b)
                    else:
                        x = 0
                        for i in range(h ** npos):
                            if (a >> i) & 1:
                                x |= b << (i * nb_new)
                        new.append(x)
            cur, npos = new, npos + k
        return Sp(cur)

    full = Sp([1 << i for i in range(m)])
    res_all = {}
    total = dict(M=0, P=0, D=0, Z=0, I0=0, Bmix=0)
    per = {}
    maxlab = 0
    hist = {}
    for j in (1, 2, 3):
        pre, post = j - 1, 3 - j
        t = tvec[j - 1]
        # factor-level spaces
        A_pre = [list(E1)] * pre                                     # full on pre positions
        if pre == 0:
            Apart = None
        # flat P (line) and B = P^perp in A on the pre positions
        if pre > 0:
            Af = Sp([flat(list(c)) for c in itertools.product(E1, repeat=pre)])
            Pf = Sp([flat(tvec[:pre])])
            Bf = perp_within(Af, Pf)
        Q = [[x] for x in tvec[j:]]                                  # lines on post positions
        def frame(kind, line=None):
            """subspaces of F_2^(h^3): kind in X_in, L0, L1, L2, L3, L5"""
            F = list(E1)
            def mk(prepart, cur):
                if pre == 0:
                    if prepart == "B":
                        return Sp()
                    parts = [cur] + Q
                else:
                    pp = {"A": ("f", Af.basis(), pre), "P": ("f", Pf.basis(), pre), "B": ("f", Bf.basis(), pre)}[prepart]
                    parts = [pp, cur] + Q
                return tens_space(parts)
            def plus(*ss):
                S = Sp()
                for s in ss:
                    for b in s.basis():
                        S.add(b)
                return S
            if kind == "X_in":
                return mk("A", [line])
            if kind == "L0":
                return mk("B", [line])
            if kind == "L1":
                return mk("B", F)
            if kind == "L2":
                return plus(mk("B", F), mk("P", [line]))
            if kind == "L3":
                return mk("A", F)
            if kind == "L5":
                tp = perp_within(Sp(E1), Sp([line])).basis()
                return plus(mk("B", F), mk("P", tp))
            raise ValueError(kind)
        cache = {}
        def fr(kind, line=None):
            k = (kind, line)
            if k not in cache:
                cache[k] = frame(kind, line)
            return cache[k]
        ZERO = Sp()
        # the x-class of triple U at this invocation: u = t_{d1} x t_{d2} x t_{d3}, varying position j-1 = U
        def uvec(U):
            tv = list(tvec)
            tv[j - 1] = ind[U]
            return flat(tv)
        S0 = a_fix[j - 1]
        nbrs = [U for U in T if nbr(S0, U)]
        d = len(nbrs)
        def lev(L, c):
            """raw level of class c ('Y' or a triple U) in a register at frame L"""
            if c == "Y":
                return L.dim()
            return L.dim() + 1 - 2 * (1 if L.reduce(uvec(c)) == 0 else 0)
        stage = dict(M=0, P=0, D=0, Bmix=0)
        def edge(name, mult, Lf, Lt, cls):
            """mult wires climb (or descend) Lf -> Lt holding the classes cls"""
            nonlocal maxlab
            if Lf.dim() <= Lt.dim():
                assert Lt.contains(Lf), name
                down = False
            else:
                assert Lf.contains(Lt), name
                down = True
            k = abs(Lt.dim() - Lf.dim())
            if k == 0:
                return
            n = len(cls)
            Dd = 0
            for c in cls:
                if c == "Y":
                    dd = k if down else 0
                else:
                    a = 1 if Lf.reduce(uvec(c)) == 0 else 0
                    b = 1 if Lt.reduce(uvec(c)) == 0 else 0
                    dd = (b - a) if not down else k - (a - b)
                assert 0 <= dd <= k, (name, c, dd, k)
                Dd += dd
                # level identity per label: lev(Lt) - lev(Lf) = (k - dd) - dd
                assert lev(Lt, c) - lev(Lf, c) == k - 2 * dd, (name, c)
            stage["M"] += mult * k
            stage["P"] += mult * k * (n - 1)
            stage["D"] += mult * Dd
            hist[n] = hist.get(n, 0) + mult * k
            maxlab = max(maxlab, n)
            per.setdefault((j, name.split()[0]), [0, 0, 0])
            per[(j, name.split()[0])][0] += mult * k
            per[(j, name.split()[0])][1] += mult * k * (n - 1)
            per[(j, name.split()[0])][2] += mult * Dd
        def gate(mult, L, before, after):
            stage["Bmix"] -= mult * (sum(lev(L, c) for c in after) - sum(lev(L, c) for c in before))
        inv = v * v                      # invocations at this stage
        S = S0
        tS = ind[S]
        X_in, Y_in = fr("X_in", tS), fr("L0", tS)
        L0S, L1, L2S, L3, L5S = fr("L0", tS), fr("L1"), fr("L2", tS), fr("L3"), fr("L5", tS)
        L4, L6, L7 = L1, L3, L3
        junk = list(nbrs)
        # ---------------- X wire S (v per invocation)
        if j == 1:      # source, holds x_S
            c = [S]
            edge("X in->2", inv * v, X_in, L2S, c); gate(inv * v, L2S, c, c)
            edge("X 2->3", inv * v, L2S, L3, c); gate(inv * v, L3, c, c)
            edge("X 3->6", inv * v, L3, L6, c); gate(inv * v, L6, c, c)
            edge("X 6->7", inv * v, L6, L7, c); gate(inv * v, L7, c, c)
        elif j == 2:    # target of the inverse invocation: x -> x - (y + x) = -y
            edge("X in->2", inv * v, X_in, L2S, [S]); gate(inv * v, L2S, [S], [S, "Y"] + junk)
            edge("X 2->3", inv * v, L2S, L3, [S, "Y"] + junk); gate(inv * v, L3, [S, "Y"] + junk, ["Y"])
            edge("X 3->6", inv * v, L3, L6, ["Y"]); gate(inv * v, L6, ["Y"], ["Y"])
            edge("X 6->7", inv * v, L6, L7, ["Y"]); gate(inv * v, L7, ["Y"], ["Y"])
        else:           # source, holds -y
            c = ["Y"]
            edge("X in->2", inv * v, X_in, L2S, c); gate(inv * v, L2S, c, c)
            edge("X 2->3", inv * v, L2S, L3, c); gate(inv * v, L3, c, c)
            edge("X 3->6", inv * v, L3, L6, c); gate(inv * v, L6, c, c)
            edge("X 6->7", inv * v, L6, L7, c); gate(inv * v, L7, c, c)
        # ---------------- Y wire S (v per invocation)
        if j == 1:      # target: y -> y + x
            c0 = ["Y"]
            edge("Y in->0", inv * v, Y_in, L0S, c0); gate(inv * v, L0S, c0, c0)
            edge("Y 0->1", inv * v, L0S, L1, c0); gate(inv * v, L1, c0, c0)
            edge("Y 1->4", inv * v, L1, L4, c0); gate(inv * v, L4, c0, ["Y", S] + junk)
            edge("Y 4->5", inv * v, L4, L5S, ["Y", S] + junk); gate(inv * v, L5S, ["Y", S] + junk, ["Y", S])
        elif j == 2:    # source, holds y + x
            c = ["Y", S]
            edge("Y in->0", inv * v, Y_in, L0S, c); gate(inv * v, L0S, c, c)
            edge("Y 0->1", inv * v, L0S, L1, c); gate(inv * v, L1, c, c)
            edge("Y 1->4", inv * v, L1, L4, c); gate(inv * v, L4, c, c)
            edge("Y 4->5", inv * v, L4, L5S, c); gate(inv * v, L5S, c, c)
        else:           # target: (y + x) + (-y) = x
            c = ["Y", S]
            edge("Y in->0", inv * v, Y_in, L0S, c); gate(inv * v, L0S, c, c)
            edge("Y 0->1", inv * v, L0S, L1, c); gate(inv * v, L1, c, c)
            edge("Y 1->4", inv * v, L1, L4, c); gate(inv * v, L4, c, c)       # junk of class Y only
            edge("Y 4->5", inv * v, L4, L5S, c); gate(inv * v, L5S, c, [S])
        # ---------------- centre wires c_k (k in [h]) and c_* (h + 1 per invocation)
        cents = [[U for U in T if kk in U] for kk in range(h)] + [list(T)]
        for members in cents:
            if j == 1:      # gathers x at L3 (row 3), read at L4 (row 4), restored at L6 (row 6)
                edge("C src->1", inv, ZERO, L1, ["Y"]); gate(inv, L1, ["Y"], ["Y"])
                edge("C 1->3", inv, L1, L3, ["Y"]); gate(inv, L3, ["Y"], ["Y"] + members)
                edge("C 3->4", inv, L3, L4, ["Y"] + members); gate(inv, L4, ["Y"] + members, ["Y"] + members)
                edge("C 4->6", inv, L4, L6, ["Y"] + members); gate(inv, L6, ["Y"] + members, ["Y"])
            elif j == 2:    # 6^-1 at L1 (with Y: += G(y + x)), 4^-1 at L3 (X reads), 3^-1 at L4 (restore), 1^-1 at L6
                edge("C src->1", inv, ZERO, L1, ["Y"]); gate(inv, L1, ["Y"], ["Y"] + members)
                edge("C 1->3", inv, L1, L3, ["Y"] + members); gate(inv, L3, ["Y"] + members, ["Y"] + members)
                edge("C 3->4", inv, L3, L4, ["Y"] + members); gate(inv, L4, ["Y"] + members, ["Y"])
                edge("C 4->6", inv, L4, L6, ["Y"]); gate(inv, L6, ["Y"], ["Y"])
            else:           # gathers -y: class Y only
                for nm, a, b in (("C src->1", ZERO, L1), ("C 1->3", L1, L3), ("C 3->4", L3, L4), ("C 4->6", L4, L6)):
                    edge(nm, inv, a, b, ["Y"]); gate(inv, b, ["Y"], ["Y"])
            edge("fin", inv, L6, full, ["Y"])
        # ---------------- side wires: physical (Y-index S, X-index U), U a neighbour of S (v d per invocation)
        for U in nbrs:
            L2U = fr("L2", ind[U])
            if j == 1:      # row 0 (with y_S, L0(S)), row 2 (+= x_U, L2(U)), row 5 (y_S reads, L5(S)), row 7 (-= x_U)
                edge("A src->0", inv * v, ZERO, L0S, ["Y"]); gate(inv * v, L0S, ["Y"], ["Y"])
                edge("A 0->2", inv * v, L0S, L2U, ["Y"]); gate(inv * v, L2U, ["Y"], ["Y", U])
                edge("A 2->5", inv * v, L2U, L5S, ["Y", U]); gate(inv * v, L5S, ["Y", U], ["Y", U])
                edge("A 5->7", inv * v, L5S, L7, ["Y", U]); gate(inv * v, L7, ["Y", U], ["Y"])
            elif j == 2:    # 7^-1 (+= y_S + x_S at L0(S)), 5^-1 (x_U reads, L2(U)), 2^-1 (-= at L5(S)), 0^-1 (L7)
                edge("A src->0", inv * v, ZERO, L0S, ["Y"]); gate(inv * v, L0S, ["Y"], ["Y", S])
                edge("A 0->2", inv * v, L0S, L2U, ["Y", S]); gate(inv * v, L2U, ["Y", S], ["Y", S])
                edge("A 2->5", inv * v, L2U, L5S, ["Y", S]); gate(inv * v, L5S, ["Y", S], ["Y"])
                edge("A 5->7", inv * v, L5S, L7, ["Y"]); gate(inv * v, L7, ["Y"], ["Y"])
            else:
                for nm, a, b in (("A src->0", ZERO, L0S), ("A 0->2", L0S, L2U), ("A 2->5", L2U, L5S), ("A 5->7", L5S, L7)):
                    edge(nm, inv * v, a, b, ["Y"]); gate(inv * v, b, ["Y"], ["Y"])
            edge("fin", inv * v, L7, full, ["Y"])
        for k in stage:
            total[k] += stage[k]
        print("h=%d stage %d: M=%d P=%d D=%d Bmix=%d" % (h, j, stage["M"], stage["P"], stage["D"], stage["Bmix"]))
    # final states: X holds -y (class Y) at full; Y holds x at U^perp; scratch class Y at full: all at level m
    N0 = v ** 3
    I0 = 3 * v * v
    W = 2 * N0 + I0 * (v * comb(h - 3, 3) + v * 3 * (h - 3) * 0 + 0)  # placeholder, recomputed below
    dn = sum(1 for U in T if nbr(T[0], U))
    W = 2 * N0 + I0 * (v * dn + h + 1)
    Mref = W * m - 2 * N0 + 2 * I0 * (h + 1) * h
    Delta = W * m - total["M"]
    ok_M = total["M"] == Mref
    ok_id = Delta == total["P"] - total["Bmix"] - 2 * total["D"]
    print("h=%d: W=%d m=%d M=%d (brute-force W m - 2N + 2L = %d, equal %s); Delta = %d" % (h, W, m, total["M"], Mref, ok_M, Delta))
    print("   raw ledger: P=%d Bmix=B=%d D=%d Z=0 I0=0 ; P - B - 2D = %d ; identity %s" %
          (total["P"], total["Bmix"], total["D"], total["P"] - total["Bmix"] - 2 * total["D"], ok_id))
    print("   per move: P/M=%.4f B/M=%.4f D/M=%.4f ; net (P-B-2D)/M = %.3e ; max labels on one move %d" %
          (total["P"] / total["M"], total["Bmix"] / total["M"], total["D"] / total["M"], Delta / total["M"], maxlab))
    print("   by stage and wire (moves, extra labels, descents):")
    for (jj, w), (M_, P_, D_) in sorted(per.items()):
        print("     stage %d %-4s M=%d P=%d D=%d P/M=%.3f" % (jj, w, M_, P_, D_, P_ / M_ if M_ else 0))
    print("   moves by number of labels carried:", dict(sorted(hist.items())))
    return dict(h=h, v=v, d=dn, m=m, W=W, M=total["M"], Delta=Delta, P=total["P"], B=total["Bmix"], D=total["D"],
                maxlab=maxlab, ok_M=ok_M, ok_identity=ok_id, hist=hist,
                per={"%d %s" % k: v_ for k, v_ in per.items()})


if __name__ == "__main__":
    out = [main(int(a)) for a in (sys.argv[1:] or ["7"])]
    os.makedirs("out", exist_ok=True)
    json.dump(out, open("out/oa_ledger_%s.json" % "_".join(sys.argv[1:] or ["7"]), "w"), indent=1, default=str)
