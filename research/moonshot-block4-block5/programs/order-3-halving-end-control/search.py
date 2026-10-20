"""search.py -- o3-build: exact template search for 2-in-7 at block 4 (5-move prefix + 2 end moves) whose first moves
use order-3 frustrated pairs (x1 + x2 against x1 - s x2, s an order-3 address map), also mixed with sign patterns.
Every prefix is tested EXACTLY with three end shapes (end2.py): PRIV (two private ends: the 640-slack branch), HALF
(+-1 halving end, E6 = E7 or E6 != E7), O3 (order-3 halving end with the template's s and s^-1, E6 = E7 or not).
Templates (generation first; all letters = progress letters of the current leads, sampled; signs +-1):
  O3S  u1 = D b1 (x1 + x2), u2 = D b2 (x1 - s x2); u3 = D(r1 u1 +- r2 u2) (both pairs advanced by ONE move);
       u4 = D(t u3 + p1 u1 +- p2 u2) (fat: depth 3 + depth 2); u5 = D(t5 u3 + t5' u4 [+ p u1/u2]) (fat)
  O3M  u1 = D b1 (x1 + x2); u2 = D(r u1 + b (x1 - s x2)) (order-3 pattern break); u3 = D(r' u2 + p u1 [+ fresh]);
       u4, u5 route sums of 1-3 earlier generators
  O3F  u1 = D b1 (x1 + x2), u2 = D b2 (x1 - s x2), u3 = D b3 (x1 - Z_v x2) or D b3 (x1 + s Z_v x2) (order-3 +
       sign-pattern gadget pieces); u4, u5 route sums with exact alignment onto a common element when possible
  O3D  two joint lineages sharing their middle move: a1 = D b1 (x1 + x2), a1' = D b1' (x1 - s x2),
       m = D(r a1 +- r' a1'), a3 = D(t m + q a1), a3' = D(t' m + q' a1')
  O3Z  Z3 orbit: u1 = D b (x1 + x2), u2 = D b (x1 - s x2), u3 = D b (x1 - s^2 x2) with ONE b, then route sums
  PL5  control: private chains 3 + 2 (must fail every end), PL6 control: 3 + 3 (PRIV must pass, 8 moves).
usage: python3 -I search.py TEMPLATE N SEED [budget_s] [out.json]"""
import sys, os, time, json, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o3lib import *
from end2 import end2, patterns, check_solution, raw_check_end2

def pick_lead(w, gen, rng):
    ld = w.leads(gen)
    return rng.choice(ld) if ld else None

def prog_letter(w, gid, rng):
    bs = w.progress(gid, want=5)
    return rng.choice(bs) if bs else None

class Builder:
    def __init__(self, w, rng):
        self.w = w; self.G = w.G; self.rng = rng
        self.gens = [w.inp(0), w.inp(1)]; self.prog = []
    def mv(self, b, cs):
        """cs: list of (coef, letter or None, j)"""
        G = self.G
        c = {}
        for (co, a, j) in cs:
            t = self.gens[j] if a is None else G.lmul(a, self.gens[j])
            c = G.add((1, c), (Fraction(co), t))
        if not c:
            return None
        u = self.w.move(b, c)
        self.gens.append(u); self.prog.append(rec_move(G, b, [(Fraction(co), a, j) for (co, a, j) in cs]))
        return len(self.gens) - 1
    def route(self, j, rng, lead=None):
        """(letter a, ) such that D1 a advances a lead term of generator j; returns a or None"""
        ld = self.w.leads(self.gens[j]) if lead is None else [lead]
        if not ld:
            return None
        (l, g) = rng.choice(ld)
        return prog_letter(self.w, g, rng)

def raw_pair(B, s, sign, a=None):
    """content x1 + sign * s x2 (s None = identity), with optional extra left letter a on x2"""
    G = B.G
    lt = s if a is None else (G.mul(a, s) if s is not None else a)
    return [(1, None, 0), (sign, lt, 1)]

def tmpl(name, w, rng, sname=None):
    G = w.G
    B = Builder(w, rng)
    keys = sorted(w.sig)
    sn = sname or rng.choice(keys); s = w.sig[sn]
    s2 = G.mul(s, s)
    sg = lambda: rng.choice([1, -1])
    first = lambda: prog_letter(w, G.ID, rng)
    if name == "PL5" or name == "PL6":
        cur = 0
        for _ in range(3):
            b = B.route(cur, rng); cur = B.mv(b, [(1, None, cur)])
        cur = 1
        for _ in range(2 if name == "PL5" else 3):
            b = B.route(cur, rng); cur = B.mv(b, [(1, None, cur)])
    elif name in ("C8O3", "C8H"):
        # control: two joint lineages T(x1 + x2), T(x1 - s x2) (C8O3) or T(x1 - x2) (C8H) with the SAME letters
        letters = []; cur = G.I
        for d in range(3):
            bb = prog_letter(w, G.intern(cur), rng); letters.append(bb); cur = G.mul(w.D1, bb, cur)
        for pair in (raw_pair(B, None, 1), raw_pair(B, s if name == "C8O3" else None, -1)):
            j = B.mv(letters[0], pair)
            for d in (1, 2):
                j = B.mv(letters[d], [(1, None, j)])
    elif name in ("G3A", "G3C", "G3S"):
        # minimal frustrated TRIPLE (tri3.py: x1 + x2, x1 + (s X_a) x2, x1 - (s X_abc) x2, ratio group of order 24, halving
        # only, no pair separates) placed at block 4 on bits (0,1,2) or (1,2,3); pieces spread over fat raw contents
        tri = rng.choice([("s012", 1, 7), ("s123", 2, 14)])
        sx = w.sig[tri[0]]
        h2 = G.mul(sx, w.X[tri[1]]); h3 = G.mul(sx, w.X[tri[2]])
        pieces = [(None, 1), (h2, 1), (h3, -1)]
        def piece(j, a):
            lt, sgn = pieces[j]
            if lt is None:
                return [(1, a, 0), (sgn, a, 1)]
            return [(1, a, 0), (sgn, (G.mul(a, lt) if a is not None else lt), 1)]
        if name == "G3A":
            # u1 carries pieces 1,2 on two frames, u2 carries pieces 2,3 (or 3,1); then the O3S fat shape
            b1 = first(); b1b = first(); b2 = first(); b2b = first()
            j1, j2, j3 = rng.sample(range(3), 3)
            u1 = B.mv(b1, piece(j1, None) + [(c, (G.mul(np.conj(b1).T, b1b, a) if a is not None else G.mul(np.conj(b1).T, b1b)), l) for (c, a, l) in piece(j2, None)])
            u2 = B.mv(b2, piece(j3, None) + [(c, (G.mul(np.conj(b2).T, b2b, a) if a is not None else G.mul(np.conj(b2).T, b2b)), l) for (c, a, l) in piece(rng.choice([j1, j2]), None)])
        elif name == "G3C":
            # three single-piece first moves would cost 3; here: u1 = piece 1, u2 = pieces 2 + 3 on two frames
            b1 = first(); b2 = first(); b2b = first()
            u1 = B.mv(b1, piece(0, None))
            u2 = B.mv(b2, piece(1, None) + [(c, (G.mul(np.conj(b2).T, b2b, a) if a is not None else G.mul(np.conj(b2).T, b2b)), l) for (c, a, l) in piece(2, None)])
        else:
            # G3S: shared base u1 = D b (x1 + x2) (piece 1), u2 = D(r u1 + pieces 2,3 fresh on two frames)
            b1 = first(); u1 = B.mv(b1, piece(0, None))
            r = B.route(u1, rng); bb = first(); bc = first()
            if r is None: return None
            ri = np.conj(r).T
            cs = [(1, None, u1)]
            for (j, bx) in ((1, bb), (2, bc)):
                cs += [(c, (G.mul(ri, bx, a) if a is not None else G.mul(ri, bx)), l) for (c, a, l) in piece(j, None)]
            u2 = B.mv(r, cs)
        r1 = B.route(u1, rng); r2 = B.route(u2, rng)
        if r1 is None or r2 is None: return None
        u3 = B.mv(r1, [(1, None, u1), (sg(), G.mul(np.conj(r1).T, r2), u2)])
        t = B.route(u3, rng); p1 = B.route(u1, rng); p2 = B.route(u2, rng)
        if any(x is None for x in (t, p1, p2)): return None
        cs = [(1, None, u3), (sg(), G.mul(np.conj(t).T, p1), u1)]
        if rng.random() < 0.6: cs.append((sg(), G.mul(np.conj(t).T, p2), u2))
        u4 = B.mv(t, cs)
        t5 = B.route(u3, rng)
        if t5 is None: return None
        d2 = [(l, g) for (l, g) in B.gens[u4] if w.depth(g) == 2 and w.depth(g) + w.dist(g) == M_]
        cs = [(1, None, u3)]
        if d2:
            a = prog_letter(w, rng.choice(d2)[1], rng)
            if a is not None: cs.append((sg(), G.mul(np.conj(t5).T, a), u4))
        if rng.random() < 0.4:
            p = B.route(u2, rng)
            if p is not None: cs.append((sg(), G.mul(np.conj(t5).T, p), u2))
        B.mv(t5, cs)
        sn = tri[0]
    elif name in ("SM", "SMO3"):
        # Lemma 2 instance (PROOF.md): a1 = D b1 c1, a1' = D b1' c2, m = D(r a1 + r' a1'), a3 = D t m, a3' = D t' m with
        # t, t' ANY progress or random letters; c1, c2 = (x1, x2) [SM] or (x1 + x2, x1 - s x2) [SMO3]
        if name == "SM":
            a1 = B.mv(first(), [(1, None, 0)]); a1p = B.mv(first(), [(1, None, 1)])
        else:
            a1 = B.mv(first(), raw_pair(B, None, 1)); a1p = B.mv(first(), raw_pair(B, s, -1))
        r = B.route(a1, rng); rp = B.route(a1p, rng)
        if any(x is None for x in (r, rp,)): return None
        m = B.mv(r, [(1, None, a1), (sg(), G.mul(np.conj(r).T, rp), a1p)])
        for _ in range(2):
            ld = w.leads(B.gens[m])
            t = prog_letter(w, rng.choice(ld)[1], rng) if (ld and rng.random() < 0.8) else w.G.mono(*mono_tuple(G, G.I)) if False else None
            if t is None:
                from fgeo import Geo
                t = Geo(4, rng.randrange(10**6)).rand_free()
            B.mv(t, [(1, None, m)])
    elif name in ("STG2", "CROSS2"):
        # skeletons outside Lemmas 2-4 (PROOF.md): raw basis = private (x1, x2) or an order-3 pair, chosen at random
        if rng.random() < 0.5:
            c1, c2 = [(1, None, 0)], [(1, None, 1)]
        else:
            c1, c2 = raw_pair(B, None, 1), raw_pair(B, s, -1)
        if name == "STG2":
            # staggered lineage + SECOND-order clean-up chain (Lemma S covers only first order at depth-3 targets)
            u1 = B.mv(first(), c1)
            r = B.route(u1, rng); bb = first()
            if r is None: return None
            ri = np.conj(r).T
            u2 = B.mv(r, [(1, None, u1)] + [(co, (G.mul(ri, bb, a) if a is not None else G.mul(ri, bb)), l) for (co, a, l) in c2])
            r3 = B.route(u2, rng)
            if r3 is None: return None
            u3 = B.mv(r3, [(1, None, u2)])
            # u4: advance u3's y2-part (depth 2 -> 3) or u2's y2-part, possibly sideways for the y1-part
            d2 = [(l, g) for (l, g) in B.gens[u3] if w.depth(g) == 2 and w.depth(g) + w.dist(g) == M_]
            t = prog_letter(w, rng.choice(d2)[1], rng) if d2 else None
            if t is None: return None
            cs = [(1, None, u3)]
            p = B.route(u2, rng)
            if p is not None and rng.random() < 0.7: cs.append((sg(), G.mul(np.conj(t).T, p), u2))
            u4 = B.mv(t, cs)
            ld = w.leads(B.gens[u4]); t5 = prog_letter(w, rng.choice(ld)[1], rng) if ld else None
            if t5 is None: return None
            cs = [(1, None, u4)]
            for j in (u3, u2):
                a = B.route(j, rng)
                if a is not None and rng.random() < 0.6: cs.append((sg(), G.mul(np.conj(t5).T, a), j))
            B.mv(t5, cs)
        else:
            # CROSS2: private first moves, a joint middle move m, then two depth-3 carriers with CROSSED private helpers
            a1 = B.mv(first(), c1); a1p = B.mv(first(), c2)
            r = B.route(a1, rng); rp = B.route(a1p, rng)
            if any(x is None for x in (r, rp,)): return None
            m = B.mv(r, [(1, None, a1), (sg(), G.mul(np.conj(r).T, rp), a1p)])
            for helper in (a1p, a1):
                t = B.route(m, rng); q = B.route(helper, rng)
                if any(x is None for x in (t, q,)): return None
                B.mv(t, [(1, None, m), (sg(), G.mul(np.conj(t).T, q), helper)])
    elif name == "CLEAN":
        # Lemma-5-directed: u1, u2 depth 1 (private or order-3 raw basis), u3 = D(r u1 +- r' u2) depth 2; the two
        # carriers u4, u5 use FRAME-SELECTIVE routes on u3: each route advances exactly one depth-2 term of u3 and no
        # other, landing on the target frames F_T1 (x1) and F_T2 (x2), so that the depth-3 parts of u4 and u5 are an
        # invertible mixing of T1 x1, T2 x2 (Lemma 5); private helpers (routes of u1 / u2) are added inside the moves
        # (outside Lemma 2); the exact end decides whether the lower-depth junk cancels.
        if rng.random() < 0.5:
            c1, c2 = [(1, None, 0)], [(1, None, 1)]
        else:
            c1, c2 = raw_pair(B, None, 1), raw_pair(B, s, -1)
        u1 = B.mv(first(), c1); u2 = B.mv(first(), c2)
        r = B.route(u1, rng); rp = B.route(u2, rng)
        if any(x is None for x in (r, rp,)): return None
        u3 = B.mv(r, [(1, None, u1), (sg(), G.mul(np.conj(r).T, rp), u2)])
        H = [(l, g) for (l, g) in B.gens[u3] if w.depth(g) == 2 and w.depth(g) + w.dist(g) == M_]
        allH = list(B.gens[u3].keys())
        def selective(rt, target):
            """route rt advances the term `target` (l, g) and no other term of u3"""
            for (l, g) in allH:
                d = w.depth(G.intern(G.mul(w.D1, rt, G.MATS[g])))
                if (l, g) == target:
                    if d != 3: return False
                elif d >= 3:
                    return False
            return True
        def find_route(l_want, frame=None, tries=4):
            cands = [(l, g) for (l, g) in H if l == l_want]
            rng.shuffle(cands)
            for (l, g) in cands:
                for rt in w.progress(g, want=8):
                    if frame is not None and G.frm(G.intern(G.mul(w.D1, rt, G.MATS[g]))) != frame:
                        continue
                    if selective(rt, (l, g)):
                        return rt, G.intern(G.mul(w.D1, rt, G.MATS[g]))
            return None
        f1 = find_route(0); f2 = find_route(1)
        if f1 is None or f2 is None: return None
        (r1, T1), (r2, T2) = f1, f2
        F1, F2 = G.frm(T1), G.frm(T2)
        f3 = find_route(0, F1); f4 = find_route(1, F2)
        if f3 is None or f4 is None: return None
        (r3, _), (r4, _) = f3, f4
        for (ra, rb, sgn_b) in ((r1, r2, 1), (r3, r4, -1)):
            cs = [(1, None, u3), (sgn_b, G.mul(np.conj(ra).T, rb), u3)]
            if rng.random() < 0.7:
                hp = B.route(rng.choice([u1, u2]), rng)
                j = rng.choice([u1, u2])
                hp = B.route(j, rng)
                if hp is not None: cs.append((sg(), G.mul(np.conj(ra).T, hp), j))
            B.mv(ra, cs)
    elif name == "O3S":
        b1 = first(); b2 = b1 if rng.random() < 0.5 else first()
        u1 = B.mv(b1, raw_pair(B, None, 1))
        u2 = B.mv(b2, raw_pair(B, s, -1))
        r1 = B.route(u1, rng); r2 = B.route(u2, rng)
        if r1 is None or r2 is None: return None
        # u3 = D r (r1' u1 +- r2' u2): one paid letter, the twist of u2 relative to u1 is free
        u3 = B.mv(r1, [(1, None, u1), (sg(), G.mul(np.conj(r1).T, r2), u2)])
        t = B.route(u3, rng); p1 = B.route(u1, rng); p2 = B.route(u2, rng)
        if any(x is None for x in (t, p1, p2,)): return None
        cs = [(1, None, u3), (sg(), G.mul(np.conj(t).T, p1), u1)]
        if rng.random() < 0.6: cs.append((sg(), G.mul(np.conj(t).T, p2), u2))
        u4 = B.mv(t, cs)
        t5 = B.route(u3, rng)
        # advance u4's depth-2 part: route a depth-2 lead of u4
        d2 = [(l, g) for (l, g) in B.gens[u4] if w.depth(g) == 2 and w.depth(g) + w.dist(g) == M_]
        if t5 is None: return None
        cs = [(1, None, u3)]
        if d2:
            a = prog_letter(w, rng.choice(d2)[1], rng)
            if a is not None: cs.append((sg(), G.mul(np.conj(t5).T, a), u4))
        if rng.random() < 0.4:
            p = B.route(u2, rng)
            if p is not None: cs.append((sg(), G.mul(np.conj(t5).T, p), u2))
        B.mv(t5, cs)
    elif name == "O3M":
        b1 = first()
        u1 = B.mv(b1, raw_pair(B, None, 1))
        r = B.route(u1, rng); bb = first()
        if r is None: return None
        u2 = B.mv(r, [(1, None, u1), (1, G.mul(np.conj(r).T, bb), 0), (-1, G.mul(np.conj(r).T, bb, s), 1)])
        r2 = B.route(u2, rng); p = B.route(u1, rng)
        if any(x is None for x in (r2, p,)): return None
        cs = [(1, None, u2), (sg(), G.mul(np.conj(r2).T, p), u1)]
        if rng.random() < 0.5:
            f2 = first(); cs += [(sg(), G.mul(np.conj(r2).T, f2), rng.choice([0, 1]))]
        u3 = B.mv(r2, cs)
        for _ in range(2):
            srcs = rng.sample(range(2, len(B.gens)), min(len(B.gens) - 2, rng.choice([1, 2, 3])))
            t = B.route(srcs[0], rng)
            if t is None: return None
            cs = [(1, None, srcs[0])]
            for j in srcs[1:]:
                a = B.route(j, rng)
                if a is not None: cs.append((sg(), G.mul(np.conj(t).T, a), j))
            B.mv(t, cs)
    elif name == "O3F":
        u1 = B.mv(first(), raw_pair(B, None, 1))
        u2 = B.mv(first(), raw_pair(B, s, -1))
        v = rng.randrange(1, N_)
        if rng.random() < 0.5:
            u3 = B.mv(first(), raw_pair(B, w.Z[v], -1))
        else:
            u3 = B.mv(first(), raw_pair(B, G.mul(s, w.Z[v]), 1))
        for _ in range(2):
            srcs = rng.sample(range(2, len(B.gens)), min(len(B.gens) - 2, rng.choice([2, 3])))
            j0 = srcs[0]
            (l0, g0) = pick_lead(w, B.gens[j0], rng)
            t = prog_letter(w, g0, rng)
            if t is None: return None
            T = G.mul(w.D1, t, G.MATS[g0])
            cs = [(1, None, j0)]
            for j in srcs[1:]:
                ld = w.leads(B.gens[j])
                done = False
                for (l, g) in ld:
                    a = w.align(T, G.MATS[g])          # land exactly on the same element T
                    if a is not None:
                        cs.append((sg(), G.mul(np.conj(t).T, a), j)); done = True; break
                if not done:
                    a = B.route(j, rng)
                    if a is not None: cs.append((sg(), G.mul(np.conj(t).T, a), j))
            B.mv(t, cs)
    elif name == "O3D":
        a1 = B.mv(first(), raw_pair(B, None, 1))
        a1p = B.mv(first(), raw_pair(B, s, -1))
        r = B.route(a1, rng); rp = B.route(a1p, rng)
        if any(x is None for x in (r, rp,)): return None
        m = B.mv(r, [(1, None, a1), (sg(), G.mul(np.conj(r).T, rp), a1p)])
        for src in (a1, a1p):
            t = B.route(m, rng); q = B.route(src, rng)
            if any(x is None for x in (t, q,)): return None
            cs = [(1, None, m), (sg(), G.mul(np.conj(t).T, q), src)]
            B.mv(t, cs)
    elif name == "O3Z":
        b = first()
        u1 = B.mv(b, raw_pair(B, None, 1))
        u2 = B.mv(b, raw_pair(B, s, -1))
        u3 = B.mv(b, raw_pair(B, s2, rng.choice([1, -1])))
        for _ in range(2):
            srcs = rng.sample(range(2, len(B.gens)), min(len(B.gens) - 2, rng.choice([2, 3])))
            t = B.route(srcs[0], rng)
            if t is None: return None
            cs = [(1, None, srcs[0])]
            for j in srcs[1:]:
                a = B.route(j, rng)
                if a is not None: cs.append((sg(), G.mul(np.conj(t).T, a), j))
            B.mv(t, cs)
    else:
        raise ValueError(name)
    return B, sn

def E_pool(w, gens, nE):
    G = w.G
    out = {}
    for gen in gens:
        for (l, gid) in gen:
            if w.dist(gid) != 1:
                continue
            E = G.mul(G.K, np.conj(G.MATS[gid]).T)
            if not (G.rowsupp(E) == 2).all():
                continue
            out.setdefault(G.key(E), E)
    vals = list(out.values())
    return vals[:nE]

def test_prefix(w, B, sn, nE=3, cap=int(os.environ.get("SEARCH_CAP", 1000))):
    G = w.G
    gens = B.gens
    res = {"PRIV": None, "HALF": [], "O3": []}
    ok, sts, sols = w.end_private(gens, cap)
    res["PRIV"] = sts
    hit = None
    if ok:
        hit = ("PRIV", None, sols)
    Es = E_pool(w, gens, nE)
    res["nE"] = len(Es)
    s = w.sig[sn]; sinv = G.mul(s, s)
    pats = [("HALF", patterns(G, s)["HALF"]), ("O3", patterns(G, s)["O3"]), ("O3", patterns(G, sinv)["O3"])]
    for i6 in range(len(Es)):
        for i7 in range(len(Es)):
            for (pn, al) in pats:
                st, sol, nc = end2(G, gens, [Es[i6], Es[i7]], al, cap)
                single = None
                if st == "no" and i6 == i7 and pn == "O3":
                    s0 = end2(G, gens, [Es[i6], Es[i7]], al, cap, exact=False, outs=(0,))[0]
                    s1 = end2(G, gens, [Es[i6], Es[i7]], al, cap, exact=False, outs=(1,))[0]
                    single = s0[0] + s1[0]
                res[pn].append((i6, i7, st, single))
                if st == "yes" and hit is None:
                    ok2 = check_solution(G, gens, [Es[i6], Es[i7]], al, sol)
                    raw = raw_check_end2(w, B.prog, gens, [Es[i6], Es[i7]], al, sol, fs=(1, 2))
                    hit = (pn, (i6, i7, ok2, {k: float('%.1e' % v) for k, v in raw.items()}), sol)
    return res, hit, Es

def main():
    name, N, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    budget = float(sys.argv[4]) if len(sys.argv) > 4 else 600
    outp = sys.argv[5] if len(sys.argv) > 5 else f"out_{name}_s{seed}.json"
    rng = random.Random(seed)
    w = W4(seed)
    t0 = time.time()
    st = {"built": 0, "fail_build": 0, "priv_one": 0, "priv_both": 0, "half_hits": 0, "o3_hits": 0,
          "single_ok_joint_no": 0, "caps": 0}
    recs = []; hits = []
    for it in range(N):
        if time.time() - t0 > budget:
            st["stopped_after"] = it; break
        if len(w.G.MATS) > 15000:
            w = W4(seed + 1000 + it)
        try:
            r = tmpl(name, w, rng)
        except Exception as e:
            r = None
        if r is None:
            st["fail_build"] += 1; continue
        B, sn = r
        st["built"] += 1
        res, hit, Es = test_prefix(w, B, sn)
        npriv = sum(1 for x in res["PRIV"] if x in ("yes", "yes0"))
        key = "|".join(sorted(res["PRIV"]))
        st.setdefault("priv_status", {}); st["priv_status"][key] = st["priv_status"].get(key, 0) + 1
        st["priv_one"] += npriv == 1; st["priv_both"] += npriv == 2
        st["half_hits"] += any(x[2] == "yes" for x in res["HALF"])
        st["o3_hits"] += any(x[2] == "yes" for x in res["O3"])
        st["single_ok_joint_no"] += any(x[3] == "yy" for x in res["HALF"] + res["O3"])
        st["caps"] += sum(1 for x in res["HALF"] + res["O3"] if x[2] == "cap")
        rec = {"it": it, "s": sn, "prog": B.prog, "PRIV": res["PRIV"], "nE": res["nE"],
               "HALF": [x[2:] for x in res["HALF"]], "O3": [x[2:] for x in res["O3"]],
               "profile": [sorted((str(k), v) for k, v in w.profile(g).items()) for g in B.gens[2:]]}
        if npriv >= 1 or any(x[3] in ("yy", "yn", "ny") for x in res["HALF"] + res["O3"]) or hit:
            recs.append(rec)
        if hit:
            pn, info, sols = hit
            rec["hit"] = [pn, str(info)]
            if pn == "PRIV":
                rec["formal"] = formal_check(w, B.gens, sols)
                rec["raw"] = raw_check(w, B.prog, B.gens, sols, fs=(1, 2))
            hits.append(rec)
            print("HIT", it, pn, info, rec.get("formal"), rec.get("raw"), flush=True)
        if it % 20 == 0:
            print(f"progress {it} {st} {time.time()-t0:.0f}s", flush=True)
    st["seconds"] = round(time.time() - t0, 1)
    print(f"search {name} N={N} seed={seed}: {st}", flush=True)
    json.dump({"args": sys.argv[1:], "stats": st, "hits": hits, "near": recs[:200]}, open(outp, "w"))

if __name__ == "__main__":
    main()
