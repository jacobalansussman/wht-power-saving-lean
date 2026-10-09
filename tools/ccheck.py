"""ccheck.py: REFERENCE CHECKER of a carrier certificate (schema carrier-cert/1), independent of the Lean files.

Reads ONLY the JSON file.  Checks the conditions of the carrier specification (a working document that is
not in this repository; the letters and section numbers below are its own):
  L1-L4 (labels: projectors recomputed from the moves), H1-H7 (shape), C (scalar identity, exact replay of all
  gates on the x-content of every role), and recomputes the price (block histograms per role kind).
scope "general" = section 10 (any ordered mix of gates x->slot, slot->slot, x->y, slot->y after the scatter);
scope "lpc"     = sections 3-4 in addition (only x->y and slot->slot in F6, H4, H7, signed supports C1-C4).

usage: /usr/bin/python3 -I -B ccheck.py <cert.json> [general|lpc]
"""
import sys, json, time


class Reject(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise Reject(msg)


def par(x):
    return bin(x).count("1") & 1


class Labels:
    """projector of every role as h row masks; a move along z needs z.z = 1 and P z = 0 (L1)."""

    def __init__(self, h, v, R, trips):
        self.h, self.v, self.R = h, v, R
        self.zero = (0,) * h
        self.P = [self.line(t) for t in trips] + [self.zero] * (R + v)
        self.moves = [0] * (2 * v + R)
        self.used = [True] * v + [False] * (R + v)
        self.blocks = {"x": {}, "y": {}, "s": {}}

    def line(self, t):
        return tuple(t if (t >> i) & 1 else 0 for i in range(self.h))

    def perp(self, t):
        return tuple((1 << i) ^ (t if (t >> i) & 1 else 0) for i in range(self.h))

    def full(self):
        return tuple(1 << i for i in range(self.h))

    def span(self, basis):
        P = list(self.zero)
        for z in basis:
            need(par(z) == 1 and not any(par(row & z) for row in P), "stored basis is not orthonormal")
            for i in range(self.h):
                if (z >> i) & 1:
                    P[i] ^= z
        return tuple(P)

    def climb(self, role, dirs, count=True):
        need(0 <= role < 2 * self.v + self.R, "role out of range")
        if not dirs:
            return
        P = list(self.P[role])
        for z in dirs:
            need(0 < z < (1 << self.h), "direction out of range")
            need(par(z) == 1, f"L1: direction {z} of role {role} is not a unit vector")
            need(not any(par(row & z) for row in P), f"L1: direction {z} of role {role} is not orthogonal to its label")
            for i in range(self.h):
                if (z >> i) & 1:
                    P[i] ^= z
        self.P[role] = tuple(P)
        self.moves[role] += len(dirs)
        if count:
            k = "x" if role < self.v else "s" if role < self.v + self.R else "y"
            self.blocks[k][len(dirs)] = self.blocks[k].get(len(dirs), 0) + 1

    def part(self, p):
        kind, role = p[0], p[1]
        need(kind in ("new", "old", "stay"), "unknown part kind")
        need(0 <= role < 2 * self.v + self.R, "role out of range")
        if kind == "stay":
            need(len(p) == 2 and self.used[role], "stay part on an unused role")
        else:
            need((kind == "new") == (not self.used[role]), f"part kind {kind} does not match the history of role {role}")
            self.used[role] = True
            self.climb(role, p[2])
        return role

    def op(self, o):
        """applies the moves of an op; L2: all parts at one label.  returns (source roles, target roles)."""
        sr = [self.part(p) for p in o["src"]]
        tg = [self.part(p) for p in o["tgt"]]
        roles = sr + tg
        need(len(set(roles)) == len(roles), "a role occurs twice in one op")
        need(len({self.P[r] for r in roles}) == 1, f"L2: gate joins unequal labels (roles {roles[:4]})")
        return sr, tg


def check(D, scope="general"):
    need(D.get("format") == "carrier-cert/1", "format")
    h, v, R = D["h"], D["v"], D["R"]
    trips = D["trips"]
    need(h >= 4 and len(trips) == v and len(set(trips)) == v, "H1: triples")
    need(all(0 < t < (1 << h) and par(t) == 1 for t in trips), "H1: a line vector is not a unit vector")
    for key, n in (("src", v), ("pieces", v), ("scat", v), ("F4", v), ("F7", v), ("ret", h), ("retLabBasis", h)):
        need(len(D[key]) == n, f"H1: length of {key}")
    need(R > 0 and len(set(D["src"])) == v and all(0 <= q < R for q in D["src"] + D["ret"]), "H1: src / ret")
    need(len(D["bank"]) == len(D["F6"]), "H5: bank table length")
    isx = lambda r: 0 <= r < v
    iss = lambda r: v <= r < v + R
    isy = lambda r: v + R <= r < 2 * v + R
    Lb = Labels(h, v, R, trips)
    con = [None] * (2 * v + R)                 # x-content: role -> {t: weight}; slots in units 1, y roles in units 1/2
    for t in range(v):
        con[t] = {t: 1}

    def add(tgt, src, num, den):
        """tgt += (num/den) src on the x-contents (exact)."""
        a = con[src]
        if not a:
            return
        f = (2 * num // den) if isy(tgt) else num
        need(isy(tgt) and den in (1, 2) or den == 1, "H7: a coefficient into a slot is not an integer")
        b = con[tgt]
        if b is None:
            b = con[tgt] = {}
        for t, w in a.items():
            need(scope != "lpc" or isy(tgt) or t not in b, "C2/C3: supports of a slot gate are not disjoint")
            x = b.get(t, 0) + f * w
            if x:
                b[t] = x
            else:
                b.pop(t, None)
    # ---- F4 (H2) --------------------------------------------------------------------------------------------
    for t, o in enumerate(D["F4"]):
        need(o["src"] == [["stay", t]] and len(o["tgt"]) == 1 and o["coef"] == [1, 1], "H2: shape of an F4 op")
        sr, tg = Lb.op(o)
        need(tg == [v + D["src"][t]], "H2: F4 target is not the slot src[t]")
        add(tg[0], t, 1, 1)
    for t in range(v):
        need(Lb.P[v + D["src"][t]] == Lb.line(trips[t]), "L3: slot src[t] is not on the line of t after F4")
    # ---- F5 (H1) --------------------------------------------------------------------------------------------
    for o in D["F5"]:
        need(len(o["src"]) == 1 and o["tgt"] and o["coef"] == [1, 1], "H1: shape of an F5 op")
        sr, tg = Lb.op(o)
        need(all(iss(r) for r in sr + tg), "H1/H2: an F5 op names a bank role")
        for r in tg:
            add(r, sr[0], 1, 1)
    # ---- the cut: retained expects, scatter -------------------------------------------------------------------
    for k in range(h):
        need(Lb.P[v + D["ret"][k]] == Lb.span(D["retLabBasis"][k]), "L3: retained slot is not at retLab after F5")
    need(all(Lb.P[v + R + S] == Lb.zero for S in range(v)), "H2: a y role has moved before the scatter")
    need(all(Lb.moves[t] == 0 for t in range(v)), "H2: an x role has moved before the scatter")
    cblocks = {}
    for bas in D["retLabBasis"]:
        cblocks[len(bas)] = cblocks.get(len(bas), 0) + 1
    for S in range(v):
        for k, num, den in D["scat"][S]:
            need(0 <= k < h and den in (1, 2), "H1: scatter entry")
            add(v + R + S, v + D["ret"][k], num, den)
    # ---- F6 (H3, H4, H5, H7) --------------------------------------------------------------------------------
    for o in D.get("F6intro", []):
        need(len(o["src"]) == 1 and not o["tgt"] and o["src"][0][2] == [], "shape of an F6intro op")
        Lb.op(o)
    Tg, Sr = set(), set()
    for o, row in zip(D["F6"], D["bank"]):
        need(len(o["src"]) == 1 and len(o["tgt"]) == 1, "H3: an F6 op is not one source and one target")
        sr, tg = Lb.op(o)
        s, t = sr[0], tg[0]
        num, den = o["coef"]
        need(row == [t, s, num, den], "H5: bank table disagrees with F6")
        need(not isy(s), "H3: an F6 gate READS a y role")
        need(not isx(t), "H3: an F6 gate goes INTO an x role")
        need(s != t and den in (1, 2) and num != 0, "H3: degenerate F6 gate")
        if scope == "lpc":
            need((isx(s) and isy(t)) or (iss(s) and iss(t)), "H3 (scope LPC): F6 gate is neither x->y nor slot->slot")
            need((den == 2 and abs(num) == 1) if isx(s) else (den == 1 and abs(num) == 1), "H7: coefficient")
            if iss(s):
                Sr.add(s)
                Tg.add(t)
        add(t, s, num, den)
    if scope == "lpc":
        need(not (Tg & Sr), "H4: a slot is both a target and a source of F6")
    # ---- F7 (H6) --------------------------------------------------------------------------------------------
    for S, o in enumerate(D["F7"]):
        pcs = D["pieces"][S]
        need(all(0 <= q < R and den == 2 and abs(num) == 1 for q, num, den in pcs), "H6/H7: piece entry")
        sr, tg = Lb.op(o)
        if pcs:
            need(sr == [v + q for q, _, _ in pcs] and tg == [v + R + S], "H6: F7 op is not (pieces of S) -> Y_S")
            need(o.get("coefs") == [[num, den] for _, num, den in pcs], "H6: F7 coefficients disagree with pieces")
        else:
            need(sr == [v + R + S] and tg == [], "H6: F7 op of a target without pieces")
        need(Lb.P[v + R + S] == Lb.perp(trips[S]), "L3: Y_S is not at t_S^perp after F7")
        for q, num, den in pcs:
            add(v + R + S, v + q, num, den)
    # ---- finals (L4) ----------------------------------------------------------------------------------------
    need([f[0] for f in D["fins"]] == list(range(2 * v + R)), "H1: final entries do not cover the roles in order")
    for r, dirs, sh, e in D["fins"]:
        if isy(r):
            need(dirs == [trips[r - v - R]], "L4: the reference move of Y_S is not along t_S")
            Lb.climb(r, dirs, count=False)
        else:
            Lb.climb(r, dirs)
        need(Lb.P[r] == Lb.full(), f"L4: role {r} does not end at the full label")
    need(all(Lb.moves[r] == (h - 1 if isx(r) else h) for r in range(2 * v + R)), "L4: move count of a role")
    need(D["N"] == R * h + v * (h - 1) + v * h, "N")
    # ---- scalar identity -------------------------------------------------------------------------------------
    for S in range(v):
        need(con[v + R + S] == {S: 2}, f"C: target {S} does not receive exactly x_S")
    bl = dict(Lb.blocks, c=cblocks)
    bl = {k: {str(r): n for r, n in sorted(d.items())} for k, d in bl.items()}
    if "blocks" in D:
        need(D["blocks"] == bl, "stored block histograms differ from the recomputed ones")
    need(sum(int(r) * n for r, n in bl["c"].items()) == D["cst"], "cst")
    return bl


def main():
    path = sys.argv[1]
    scope = sys.argv[2] if len(sys.argv) > 2 else "lpc"
    t0 = time.time()
    D = json.load(open(path))
    try:
        bl = check(D, scope)
    except Reject as e:
        print(f"REJECTED ({scope}): {e}")
        sys.exit(1)
    h, v, R = D["h"], D["v"], D["R"]
    tot = {k: sum(int(r) * n for r, n in d.items()) for k, d in bl.items()}
    nb = {k: sum(d.values()) for k, d in bl.items()}
    print(f"ACCEPTED ({scope}): {path}\n  h={h} v={v} R={R} ops F4/F5/F6/F7={len(D['F4'])}/{len(D['F5'])}/"
          f"{len(D['F6'])}/{len(D['F7'])}; moves x/y/s/c={tot} (expected {v * (h - 1)}, {v * (h - 1)}, {R * h}, "
          f"{(h - 1) ** 2 + h}); blocks x/y/s/c={nb}\n  x blocks {bl['x']}\n  y blocks {bl['y']}\n  slot blocks "
          f"{bl['s']}\n  copy blocks {bl['c']}  [{time.time() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
