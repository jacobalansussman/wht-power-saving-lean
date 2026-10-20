"""run_fam.py -- o3-shapeR: shape-R prefix families, each tested with the exact eps-quotient test (shapeR.py) against
EVERY target pair (F1, F2) (depth 3, one move from K, adjacent to a depth-2 term of u3) and ALL u4, u5 at once.
Families (u1 = D b1 y1, u2 = D b2 y2 with geodesic b_l; h1 = D r D b1, h2 = D r' D b2 geodesic, depth 2):
  A    Theorem A control: u3 = D(r u1 + s r' u2)                                  (must never pass)
  RAWS raw extras placed on the depth-1 junk of selective letters t, t' (o3-close's R0/R1, probeR)
  RAWR raw extras w1 y1, w2 y2 on random classes (1-2 per input), random signs
  SECS second letter on u2 (u1) placing a u3 frame on the depth-2 junk of t (t')
  SECR random second letters on u1 and/or u2 (any class), random signs
  COS  second letter in the SAME class as r (r2 = h r, h in H_D): two terms on one frame (zero-divisor coefficients)
  MIX  random mixture: 1-3 letters on u1, 1-3 on u2, 0-2 raw per input, coefficients in {+-1, +-2, +-1/2}
usage: python3 -I run_fam.py FAMILY NPREFIX SEED OUT.jsonl"""
import sys, os, time, json, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shapeR import *

fam = sys.argv[1]; NP = int(sys.argv[2]); seed = int(sys.argv[3]); outp = sys.argv[4]
W = World(seed); G = W.G; rng = random.Random(seed * 7 + 1)
y1, y2 = gen_inp(0), gen_inp(1)
HD = [b for b in (W.rand_letter() for _ in range(4000)) if fk(G.mul(W.D, b)) == fk(W.D)]

def geo_letter(M, want, tries=3000):
    for _ in range(tries):
        b = W.rand_letter()
        if W.st(G.mul(W.D, b, M)) == want:
            return b
    return None

def cls_of_frame(M, fkey):
    """a letter b (class rep times H_D element) with frame(D b M) == fk"""
    for b in W.classes:
        if fk(G.mul(W.D, b, M)) == fkey:
            return b
    return None

def sel_letters(h1, h2):
    """classes t: D t h1 at (3,1) and depth(D t h2) <= 2; returns list of (class index, junk matrix)"""
    out = []
    for N, b in enumerate(W.classes):
        if W.st(G.mul(W.D, b, h1)) == (3, 1):
            J = G.mul(W.D, b, h2)
            if W.st(J)[0] <= 2:
                out.append((N, J))
    return out

def sgn():
    return rng.choice((1, -1))

def coef():
    return rng.choice((1, -1, 2, -2, Fr(1, 2), Fr(-1, 2)))

def build():
    b1 = geo_letter(G.I, (1, 3)); b2 = geo_letter(G.I, (1, 3))
    if any(x is None for x in (b1, b2)):
        return None
    a1, a2 = G.mul(W.D, b1), G.mul(W.D, b2)
    r = geo_letter(a1, (2, 2)); rp = geo_letter(a2, (2, 2))
    if any(x is None for x in (r, rp)):
        return None
    h1, h2 = G.mul(W.D, r, a1), G.mul(W.D, rp, a2)
    u1 = move(W, b1, y1); u2 = move(W, b2, y2)
    parts = [(1, lmul(G, r, u1, G.I)), (sgn(), lmul(G, rp, u2, G.I))]
    info = {}
    if fam == "A":
        pass
    elif fam in ("RAWS", "SECS"):
        S1 = sel_letters(h1, h2); S2 = sel_letters(h2, h1)
        if not S1 or not S2:
            return "nosel"
        (N1, J2) = rng.choice(S1); (N2, J1) = rng.choice(S2)
        d2, d1 = W.st(J2)[0], W.st(J1)[0]
        info["junkdepth"] = [d1, d2]
        for (J, l, a, d) in ((J1, 0, a1, d1), (J2, 1, a2, d2)):
            fkj = fk(J)
            yl = y1 if l == 0 else y2; ul = u1 if l == 0 else u2
            if fam == "RAWS" and d == 1:
                w = cls_of_frame(G.I, fkj)
                parts.append((sgn(), lmul(G, w, yl, G.I)))
            elif fam == "SECS" and d == 2:
                w = cls_of_frame(a, fkj)
                if w is None:
                    info.setdefault("nonadj", 0); info["nonadj"] += 1
                    continue
                parts.append((sgn(), lmul(G, w, ul, G.I)))
        if len(parts) == 2:
            return "noextra"
    elif fam == "RAWR":
        for (l, yl) in ((0, y1), (1, y2)):
            for _ in range(rng.choice((1, 1, 2))):
                parts.append((sgn(), lmul(G, rng.choice(W.classes), yl, G.I)))
    elif fam == "SECR":
        for ul in rng.choice(([u1], [u2], [u1, u2])):
            parts.append((sgn(), lmul(G, rng.choice(W.classes), ul, G.I)))
    elif fam == "COS":
        which = rng.choice((0, 1, 2))
        if which in (0, 2):
            parts.append((sgn(), lmul(G, G.mul(rng.choice(HD), r), u1, G.I)))
        if which in (1, 2):
            parts.append((sgn(), lmul(G, G.mul(rng.choice(HD), rp), u2, G.I)))
    elif fam == "MIX":
        parts = [(coef(), lmul(G, r, u1, G.I)), (coef(), lmul(G, rp, u2, G.I))]
        for _ in range(rng.randrange(0, 3)):
            parts.append((coef(), lmul(G, rng.choice(W.classes + [G.mul(h, r) for h in HD[:3]]), u1, G.I)))
        for _ in range(rng.randrange(0, 3)):
            parts.append((coef(), lmul(G, rng.choice(W.classes + [G.mul(h, rp) for h in HD[:3]]), u2, G.I)))
        for (yl) in (y1, y2):
            for _ in range(rng.randrange(0, 3)):
                parts.append((coef(), lmul(G, rng.choice(W.classes), yl, G.I)))
    elif fam in ("MIXED", "MIXCOS"):
        # the case Lemmas R1/R2 leave open: on each input >= 2 depth-2 frames AND >= 1 raw frame
        D2 = [b for b in W.classes]
        for (l, ul, yl, a, rr) in ((0, u1, y1, a1, r), (1, u2, y2, a2, rp)):
            ks = [b for b in W.classes if W.st(G.mul(W.D, b, a))[0] == 2 and fk(G.mul(W.D, b, a)) != fk(G.mul(W.D, rr, a))]
            for b2 in rng.sample(ks, rng.choice((1, 2))):
                parts.append((coef(), lmul(G, b2, ul, G.I)))
            for _ in range(rng.choice((1, 2))):
                parts.append((coef(), lmul(G, rng.choice(W.classes), yl, G.I)))
            if fam == "MIXCOS":
                parts.append((coef(), lmul(G, G.mul(rng.choice(HD), rr), ul, G.I)))
    elif fam in ("TYPEM", "TYPEM1", "TYPEMR"):
        # the case no lemma covers: a raw frame W adjacent to BOTH depth-2 frames of u3 that it can touch, i.e. W
        # adjacent to H = frame(h_l) and H' = the other common neighbour of N_l and W, with a letter on u_l giving H'.
        # TYPEM: both inputs; TYPEM1: input 1 only (input 2 plain).
        for (l, ul, yl, a, hh) in ((0, u1, y1, a1, h1), (1, u2, y2, a2, h2)):
            if fam in ("TYPEM1", "TYPEMR") and l == 1:
                if fam == "TYPEMR":
                    parts.append((coef(), lmul(G, rng.choice(W.classes), yl, G.I)))
                continue
            nbH = W.nbrs(hh)[1]
            Ws = [X for X in nbH if abs(W.dep[X] - 1) < 1e-9 and X != fk(a)]
            if not Ws:
                return "noW"
            X = rng.choice(Ws)
            w = cls_of_frame(G.I, X)
            nbN = set(W.nbrs(a)[1]); nbX = None
            Dw = G.mul(W.D, w); nbX = set(W.nbrs(Dw)[1])
            Hs = [f for f in nbN & nbX if abs(W.dep[f] - 2) < 1e-9 and f != fk(hh)]
            if not Hs:
                return "noH2"
            b2 = cls_of_frame(a, Hs[0])
            parts.append((coef(), lmul(G, G.mul(rng.choice(HD), b2), ul, G.I)))
            parts.append((coef(), lmul(G, G.mul(rng.choice(HD), w), yl, G.I)))
    u3 = move(W, G.I, comb(*parts))
    return [y1, y2, u1, u2, u3], info

stats = collections.Counter(); t0 = time.time(); fo = open(outp, "w"); built = 0; it = 0
while built < NP and it < 20 * NP:
    it += 1
    b = build()
    if b is None or isinstance(b, str):
        stats[b or "nobuild"] += 1; continue
    gens, info = b
    u3 = gens[4]
    res, n1, n2 = test_prefix(W, gens)
    built += 1
    # relaxed ladder: lowest level with a passing pair (3: junk at depth < 3 ignored; 2: depth < 2 ignored; 0 exact)
    best = None
    if any(x["pass"] for x in res):
        best = 0
    else:
        r3, _, _ = test_prefix(W, gens, lev=3)
        if any(x["pass"] for x in r3):
            best = 3
            r2, _, _ = test_prefix(W, gens, lev=2)
            if any(x["pass"] for x in r2):
                best = 2
    stats[("bestlev", best)] += 1
    npass = sum(x["pass"] for x in res)
    d4max = max((x["d4"][0] for x in res), default=-1); d5max = max((x["d5"][0] for x in res), default=-1)
    stats["built"] += 1; stats["pass_prefix"] += npass > 0; stats[("d4max", d4max)] += 1; stats[("d5max", d5max)] += 1
    stats[("nT", n1, n2)] += 1
    rec = {"fam": fam, "seed": seed, "k": built, "u3terms": len(u3), "nT": [n1, n2], "npairs": len(res), "npass": npass,
           "d4": [x["d4"] for x in res], "d5": [x["d5"] for x in res], "info": info, "bestlev": best}
    fo.write(json.dumps(rec) + "\n"); fo.flush()
    if len(G.MATS) > 20000:
        G.MATS.clear(); G.IDX.clear(); G.FREE.clear(); G.FRM.clear(); G.KID = G.intern(G.K); G.ID = G.intern(G.I)
    if len(W.cache) > 200000:
        W.cache.clear()
print(fam, "seed", seed, "built", built, "iters", it, "time", round(time.time() - t0, 1), "s")
for k, v in sorted(stats.items(), key=lambda kv: str(kv[0])):
    print("  ", k, v)
