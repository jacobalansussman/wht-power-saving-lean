#!/usr/bin/env python3
"""p-wgen (round wht45/lever, 19 Oct 2026): a gcert/1w word (p-lower FORMAT.md: #3945's BORROW words with warps,
twisted adds, zero buffers, R = 408) -> the Lean DATA and KERNEL-CHECK modules for p-warpchk's warp checker
(`Work/KLowW/{Data,LabK,ScalK}.lean`, namespace KLW: events `KLW.WEv`, warp table `KLW.WOp`, label replay
`KLW.segW`, free-warp-algebra scalar replay `KLW.scalW`), after an exact PYTHON MIRROR of both kernel functions.

The mirror is own code, written from the Lean text of p-warpchk's LabK.lean / ScalK.lean (and GLab.gateK of
Work/GCert/Labels/Def.lean); the reference semantics is p-lower's lower_w.check_frames / full_replay (FORMAT.md 5),
which the mirror is compared against (frames at the end, hub endpoint).  Nothing is written if the mirror rejects.

Program (gcert/1w `order`): D0, A[0:nP] (P), C0, A[nP:], scatter, B, U.  Events: `out` -> `.g (.out ..)`, `warp` ->
`.w r op f0 f1`, `zalloc z` -> `.za z 1`, `zfree z` -> `.zf z`.  Registers 0..631 and the 91 zero buffers 632..722.

usage: wgen.py <word.gcert1w.json.gz> <NAME> <out dir> [--seg N events per segment (default 2000)]
       [--mirror-only] [--mod-bytes B data bytes per module (default 700000)]
Writes <out>/Work/GCert/Data/Gen/<NAME>/D<i>.lean (data), <out>/Work/GCert/Data/Gen/<NAME><Check>.lean (one kernel
evaluation each), <out>/Work/GCert/Data/Gen/<NAME>Hist.lean, <NAME>Fp.lean, and <out>/logs/wgen-<NAME>.json."""
import sys, os, json, gzip, hashlib, argparse
from fractions import Fraction as Q
from math import lcm
from collections import Counter


class Reject(Exception):
    pass


def need(c, msg):
    if not c:
        raise Reject(msg)


# ------------------------------------------------------------------ GF(2) and frame codes
def colK(cols, m):
    """KLW.colK: XOR of the columns of the bits of m (bit i -> cols[i]; bits beyond the list are ignored)"""
    r, i = 0, 0
    while m and i < len(cols):
        if m & 1:
            r ^= cols[i]
        m >>= 1
        i += 1
    return r


def inv_cols(cols, h):
    """the columns of G^-1 (G given by its columns), by Gauss-Jordan over GF(2); None if singular"""
    rows = [(cols[i], 1 << i) for i in range(h)]       # pairs (G e_i, e_i): reduce the left to unit vectors
    basis = {}
    for a, b in rows:
        for p in sorted(basis, reverse=True):
            if a >> p & 1:
                a ^= basis[p][0]
                b ^= basis[p][1]
        if not a:
            return None
        p = a.bit_length() - 1
        for q in list(basis):
            if basis[q][0] >> p & 1:
                basis[q] = (basis[q][0] ^ a, basis[q][1] ^ b)
        basis[p] = (a, b)
    return [basis[i][1] for i in range(h)]              # e_i = G(b): Gi(e_i) = b


def code(h, basis):
    """GLab code: dim + 2^8 M + 2^(8+h) P (M pivot mask, the vector of pivot p at digit p)"""
    M = P = 0
    for b in basis:
        p = b.bit_length() - 1
        M |= 1 << p
        P |= b << (h * p)
    return len(basis) + (M << 8) + (P << (8 + h))


def in_span(basis, x):
    """x in the span of a reduced echelon basis (pivots = leading bits)"""
    for b in sorted(basis, reverse=True):
        if x >> (b.bit_length() - 1) & 1:
            x ^= b
    return x == 0


def maps_into(cols, hm, fa, fb):
    """KLW.trK half: the image of every basis vector of fa under cols (masked) lies in the span of fb"""
    return all(in_span(fb, colK(cols, u) & hm) for u in fa)


# ------------------------------------------------------------------ events
def co(num, den):
    need(den > 0, "denominator")
    return (num < 0, abs(num), den)


def events(c):
    """the program as event lists per phase (tuples): ('g', ('out', f, s, [(t, co)], ex)) | ('w', r, op, f0, f1)
    | ('za', z, 1) | ('zf', z)"""
    def ev(g):
        k = g[0]
        if k == "out":
            return ("g", ("out", g[1], g[2], [(t, co(a, b)) for t, a, b in g[3]], list(g[4])))
        if k == "in":
            return ("g", ("in", g[1], g[2], [(s, co(a, b)) for s, a, b in g[3]]))
        if k == "warp":
            return ("w", g[1], g[2], g[3], g[4])
        if k == "zalloc":
            return ("za", g[1], 1)
        if k == "zfree":
            return ("zf", g[1])
        raise Reject("unknown gate kind %r" % (k,))
    need(c["order"] == ["D0", "A[:nP]", "C0", "A[nP:]", "scat", "B", "U"], "unexpected program order")
    nP = c["nP"]
    ph = [("D0", c["D0"]), ("P", c["A"][:nP]), ("C0", c["C0"]), ("A", c["A"][nP:]), ("B", c["B"]), ("U", c["U"])]
    return [(nm, [ev(g) for g in gs]) for nm, gs in ph]


def gate_regs(g):
    if g[0] == "out":
        return [g[2]] + [t for t, _ in g[3]] + list(g[4])
    return [g[2]] + [s for s, _ in g[3]]


def gate_adds(g):
    """GXD.Gate.adds: (tgt, src, co) in order"""
    if g[0] == "out":
        return [(t, g[2], k) for t, k in g[3]]
    return [(g[2], s, k) for s, k in g[3]]


# ------------------------------------------------------------------ label mirror (KLW.evK / segW, GLab.gateK)
class Lab:
    def __init__(s, c, n):
        s.h, s.n = c["h"], n
        s.hm = (1 << s.h) - 1
        s.fr = [list(f) for f in c["frames"]]
        s.dim = [len(f) for f in s.fr]
        s.wt = [w["G_columns"] for w in c["warps"]]
        need([w["op"] for w in c["warps"]] == list(range(1, len(c["warps"]) + 1)), "warp ops not 1..")
        s.wti = []
        for i, G in enumerate(s.wt):
            Gi = inv_cols(G, s.h)
            need(Gi is not None, "warp %d singular" % (i + 1))
            need(all(colK(G, Gi[j]) & s.hm == 1 << j and colK(Gi, G[j]) & s.hm == 1 << j for j in range(s.h)),
                 "wtK fails at op %d" % (i + 1))
            s.wti.append(Gi)
        s.most = 1

    def gate(s, g, S, rs):
        f = g[1]
        need(f < len(s.fr), "frame id out of the table")
        db = s.dim[f]
        A = []
        climbers = 0
        for r in gate_regs(g):
            need(r < s.n, "register out of range")
            a = S[r]
            if a != f:
                da = s.dim[a]
                need(da < db, ("climb not to a larger dimension", g))
                rs.append(db - da)
                A.extend(s.fr[a])
                S[r] = f
                climbers += 1
        s.most = max(s.most, climbers)
        need(all(in_span(s.fr[f], u) for u in A), ("climb not nested", g))

    def event(s, e, S, rs):
        k = e[0]
        if k == "g":
            s.gate(e[1], S, rs)
        elif k == "w":
            _, r, a, f0, f1 = e
            need(r < s.n and S[r] == f0 and a != 0 and abs(a) <= len(s.wt), ("warp event", e, S[r] if r < s.n else None))
            G, Gi = s.wt[abs(a) - 1], s.wti[abs(a) - 1]
            ca, cb = (f0, f1) if a > 0 else (f1, f0)
            need(maps_into(G, s.hm, s.fr[ca], s.fr[cb]) and maps_into(Gi, s.hm, s.fr[cb], s.fr[ca]), ("transK", e))
            S[r] = f1
        elif k == "za":
            need(e[1] < s.n, "za register")
            S[e[1]] = e[2]
        elif k == "zf":
            need(e[1] < s.n, "zf register")
        else:
            raise Reject(k)


# ------------------------------------------------------------------ scalar mirror (KLW.evS / scalW)
def enc(a):
    return 2 * a if a > 0 else 2 * (-a - 1) + 1


def encN(a):
    return enc(-a)


def wcode(B, a, cw):
    return cw // B if cw % B == encN(a) else enc(a) + B * cw


class Scal:
    """rows as dicts key -> Fraction value (exact); the Lean data is z = value * D (asserted integral)"""

    def __init__(s, nr, N, B):
        s.nr, s.N, s.B = nr, N, B
        s.dens = 1
        s.maxenc = 0

    def scaled(s, row, ng, num, den):
        out = {}
        for k, x in row.items():
            y = x * num / den
            if ng:
                y = -y
            s.dens = lcm(s.dens, y.denominator)
            out[k] = y
        return out

    def add(s, R, t, src, cf):
        ng, num, den = cf
        need(t < s.nr and src < s.nr and t != src and den != 0, ("addS", t, src))
        rs = s.scaled(R[src], ng, num, den)
        rt = R[t]
        for k, y in rs.items():
            z = rt.get(k, 0) + y
            if z:
                rt[k] = z
                s.dens = lcm(s.dens, z.denominator)
            else:
                rt.pop(k, None)

    def event(s, e, R, o):
        k = e[0]
        if k == "g":
            for t, src, cf in gate_adds(e[1]):
                if not o:
                    s.add(R, t, src, cf)
                else:
                    s.add(R, src, t, (not cf[0], cf[1], cf[2]))
        elif k == "w":
            _, r, a, _, _ = e
            need(r < s.nr and a != 0, "warp")
            s.maxenc = max(s.maxenc, enc(a), encN(a))
            need(enc(a) < s.B and encN(a) < s.B, "B too small")
            out = {}
            for key, x in R[r].items():
                src, cw = key % s.N, key // s.N
                out[src + s.N * wcode(s.B, a, cw)] = x
            R[r] = out
        elif k in ("za", "zf"):
            need(e[1] < s.nr, k)
            R[e[1]] = {}
        else:
            raise Reject(k)


def word_of(key, N, B):
    """(source, word) of a key (word as a tuple of signed letters, first = leading letter)"""
    src, cw = key % N, key // N
    w = []
    while cw:
        d = cw % B
        w.append(d // 2 if d % 2 == 0 else -(d // 2 + 1))
        cw //= B
    return src, tuple(w)


# ------------------------------------------------------------------ Lean text
def ll(xs):
    return "[" + ", ".join(xs) + "]"


def lint(z):
    return str(z) if z >= 0 else "Int.negSucc %d" % (-z - 1)


class Emit:
    def __init__(s):
        s.defs = []
        s.co = {}

    def cf(s, k):
        if k not in s.co:
            s.co[k] = "k%d" % len(s.co)
        return s.co[k]

    def add(s, nm, typ, body):
        s.defs.append("noncomputable def %s : %s := %s" % (nm, typ, body))

    def lst(s, nm, typ, items, fmt, per=250):
        parts = []
        for i in range(0, max(len(items), 1), per):
            pn = "%s_%d" % (nm, i // per)
            parts.append(pn)
            s.add(pn, "List (%s)" % typ, ll(map(fmt, items[i:i + per])))
        s.add(nm, "List (%s)" % typ, " ++ (".join(parts) + ")" * (len(parts) - 1))

    def gate(s, g):
        pr = lambda x: "(%d, %s)" % (x[0], s.cf(x[1]))
        if g[0] == "out":
            return ".out %d %d %s %s" % (g[1], g[2], ll(map(pr, g[3])), ll(map(str, g[4])))
        return ".inn %d %d %s" % (g[1], g[2], ll(map(pr, g[3])))

    def ev(s, e):
        k = e[0]
        if k == "g":
            return ".g (%s)" % s.gate(e[1])
        if k == "w":
            return ".w %d (%s) %d %d" % (e[1], lint(e[2]), e[3], e[4])
        if k == "za":
            return ".za %d %d" % (e[1], e[2])
        return ".zf %d" % e[1]

    def trie(s, nm, vals, d, low=8):
        """SSC.Trie of depth d (key bit 0 first: left = even keys), literals of 2^low leaves"""
        vals = list(vals) + [0] * ((1 << d) - len(vals))
        cnt = [0]

        def lit(vs):
            return ".leaf %d" % vs[0] if len(vs) == 1 else ".node (%s) (%s)" % (lit(vs[0::2]), lit(vs[1::2]))

        def go(vs, depth, path):
            if depth <= low:
                pn = "%s_%d" % (nm, cnt[0])
                cnt[0] += 1
                s.add(pn, "SSC.Trie", lit(vs))
                return pn
            l, r = go(vs[0::2], depth - 1, path + "0"), go(vs[1::2], depth - 1, path + "1")
            pn = nm if not path else "%s_n%s" % (nm, path)
            s.add(pn, "SSC.Trie", ".node %s %s" % (l, r))
            return pn
        root = go(vals, d, "")
        if root != nm:
            s.add(nm, "SSC.Trie", root)

    def rtrie(s, nm, rows, d, low=5):
        """KLW.RT of depth d: leaf i = rows[i] (list of (key, z)), [] beyond"""
        rows = list(rows) + [[]] * ((1 << d) - len(rows))
        cnt = [0]
        row = lambda r: ll("(%d, %s)" % (k, lint(z)) for k, z in r)

        def lit(rs):
            return ".leaf %s" % row(rs[0]) if len(rs) == 1 else ".node (%s) (%s)" % (lit(rs[0::2]), lit(rs[1::2]))

        def go(rs, depth, path):
            if depth <= low:
                pn = "%s_%d" % (nm, cnt[0])
                cnt[0] += 1
                s.add(pn, "KLW.RT", lit(rs))
                return pn
            l, r = go(rs[0::2], depth - 1, path + "0"), go(rs[1::2], depth - 1, path + "1")
            pn = nm if not path else "%s_n%s" % (nm, path)
            s.add(pn, "KLW.RT", ".node %s %s" % (l, r))
            return pn
        root = go(rows, d, "")
        if root != nm:
            s.add(nm, "KLW.RT", root)


def bits(n):
    d = 0
    while (1 << d) < n:
        d += 1
    return d


# ------------------------------------------------------------------ fingerprint (mirror of PWG.fp, Work/PWGen/Fp.lean)
FNV_OFF, FNV_P, M64 = 14695981039346656037, 1099511628211, 1 << 64


def mix(hh, x):
    return ((hh ^ x) * FNV_P) % M64


def fp_nats(c, evs, wt, wti):
    """the Nat stream of PWG.fpWord: segments' events, then wt (G, Gi), ret, scat rows, start, final"""
    out = []
    def co_(k):
        out.extend([1 if k[0] else 0, k[1], k[2]])
    def lst(xs, f):
        out.append(len(xs))
        for x in xs:
            f(x)
    def gate(g):
        if g[0] == "out":
            out.extend([0, g[1], g[2]])
            lst(g[3], lambda p: (out.append(p[0]), co_(p[1])))
            lst(g[4], out.append)
        else:
            out.extend([1, g[1], g[2]])
            lst(g[3], lambda p: (out.append(p[0]), co_(p[1])))
    def ev(e):
        if e[0] == "g":
            out.append(0); gate(e[1])
        elif e[0] == "w":
            out.extend([2, e[1], enc(e[2]), e[3], e[4]])
        elif e[0] == "za":
            out.extend([3, e[1], e[2]])
        else:
            out.extend([4, e[1]])
    lst(evs, lambda seg: lst(seg, ev))
    lst(list(zip(wt, wti)), lambda p: (lst(p[0], out.append), lst(p[1], out.append)))
    lst(sorted(c["ret"]), lambda r: out.extend([r[1], r[2]]))
    lst(c["scat"]["table"], lambda row: lst(row, lambda x: (out.append(x[0]), co_(co(x[1], x[2])))))
    lst(c["start"], out.append)
    lst(c["final"], out.append)
    return out


def fnv(ns):
    hh = FNV_OFF
    for x in ns:
        hh = mix(hh, x)
    return hh


# ------------------------------------------------------------------ main
HEAD = """/-! gcert/1w word `%(name)s` (#2979 provider %(prov)s; h = 9, v = 112, R = 408, 91 zero buffers): %(what)s.
Source: %(src)s (sha256 %(sha)s).  Generated by p-wgen's tools/wgen.py, x-words copy for #2979 (round wht45/lever2979; do not edit). -/
set_option maxRecDepth 100000
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cert")
    ap.add_argument("name")
    ap.add_argument("out")
    ap.add_argument("--seg", type=int, default=2000)
    ap.add_argument("--mod-bytes", type=int, default=700000)
    ap.add_argument("--mirror-only", action="store_true")
    ap.add_argument("--expect", default='{"1": 2482, "2": 712, "3": 218, "4": 109, "5": 20, "8": 8}',
                    help="expected local profile (default: #3945's)")
    a = ap.parse_args()
    raw = open(a.cert, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    c = json.loads(gzip.decompress(raw))
    need(c.get("format") == "gcert/1w" and c["copies"] == 0, "not a gcert/1w file without copies")
    NAME = a.name
    h, v, R, nZ = c["h"], c["v"], c["R"], c["nZ"]
    n0 = 2 * v + R
    n = n0 + nZ
    need([z for z, _ in c["zero_buffers"]] == list(range(n0, n)), "zero buffer numbering")
    phases = events(c)
    # ---- segments: cut every `seg` events and at every phase boundary; the scatter sits between P..A and B..U
    segs, segph = [], []
    for nm, evs in phases:
        for i in range(0, len(evs), a.seg):
            segs.append(evs[i:i + a.seg])
            segph.append(nm)
    nA = sum(1 for p in segph if p in ("D0", "P", "C0", "A"))
    # ---- label mirror
    L = Lab(c, n)
    S = list(c["start"]) + [0] * nZ
    states, rss = [list(S)], []
    for i, sg in enumerate(segs):
        rs = []
        for e in sg:
            L.event(e, S, rs)
        rss.append(rs)
        states.append(list(S))
    sA = states[nA]
    levels = {int(t): f for t, f in c["levels"].items()}
    need(all(sA[v + t] == levels.get(t, 0) for t in range(v)), "y not at its level at the scatter")
    need(all(sA[s] == f for _, s, f in c["ret"]), "a retained total is not at its frame at the scatter")
    cen = [L.dim[f] for _, _, f in c["ret"]]
    need(states[-1][:n0] == c["final"], "the label replay does not end at the final frames")
    hist = Counter(sum(rss, [])) + Counter(cen)
    need({str(r): x for r, x in sorted(hist.items())} == c["hub_local_profile"], ("histogram", dict(hist)))
    need(dict(sorted(hist.items())) == {int(k): x for k, x in json.loads(a.expect).items()}, ("not the expected local profile", a.expect))
    # ---- scalar mirror, both orientations, exact
    N, B = n, 256
    scal = {}
    for o in (False, True):
        Sc = Scal(n, N, B)
        Rw = [{r: Q(1)} for r in range(n)]          # identity rows for ALL registers (zero buffers too: `.za` erases them)
        snaps = [[dict(x) for x in Rw]]
        tot = None
        for i, sg in enumerate(segs):
            if i == nA:
                # the scatter: y_t += coef * total_k (totals read at the cut; never written by the scatter)
                ret = {k: s for k, s, _ in c["ret"]}
                sc_ev = [("g", ("in", 0, v + t, [(ret[k], co(p, q))])) for t in range(v) for k, p, q in c["scat"]["table"][t]]
                for e in sc_ev:
                    Sc.event(e, Rw, o)
                snaps.append([dict(x) for x in Rw])
            for e in sg:
                Sc.event(e, Rw, o)
            snaps.append([dict(x) for x in Rw])
        # end state, as FORMAT.md 5.3 (forward) and its inverse transpose (dual)
        ue = {r: Q(p, q) for r, p, q in c["units_end"]}
        bad = []
        for r in range(n0):
            if not o:
                want = {r: 1 / ue.get(r, Q(1))}
                if v <= r < 2 * v:
                    want[r - v] = Q(1)
            else:
                want = {r: ue.get(r, Q(1))}
                if r < v:
                    want[v + r] = Q(-1)
            if Rw[r] != want:
                bad.append(r)
        need(not bad, ("endpoint fails", o, bad[:8]))
        need(all(not Rw[z] for z in range(n0, n)), "zero buffers not empty at the end")
        D = Sc.dens
        scal[o] = dict(snaps=snaps, D=D, maxenc=Sc.maxenc, sc_ev=sc_ev)
    D = lcm(scal[False]["D"], scal[True]["D"])
    # integrality of every stored and every scaled coefficient times D is implied by D = lcm of their denominators
    info = dict(name=NAME, cert=a.cert, cert_sha256=sha, provider=c["provider"], n=n, segments=len(segs), nA=nA,
                seg_phases=segph, seg_lengths=[len(x) for x in segs], ranks=dict(sorted(hist.items())),
                D=D, B=B, N=N, K=h * L.most, most=L.most, maxenc=max(scal[False]["maxenc"], scal[True]["maxenc"]),
                max_row_entries=[max(sum(len(x) for x in sn) for sn in scal[o]["snaps"]) for o in (False, True)])
    # ---- fingerprint
    ns = fp_nats(c, segs, L.wt, L.wti)
    info["fp"] = fnv(ns)
    info["fp_stream_length"] = len(ns)
    print(json.dumps(info, indent=1))
    if a.mirror_only:
        return
    # ---- Lean text
    E = Emit()
    tag = dict(name=NAME, prov=c["provider"], src=os.path.basename(a.cert), sha=sha)
    for i, sg in enumerate(segs):
        E.lst("ev%d" % i, "KLW.WEv", sg, E.ev, per=150)
        E.lst("rs%d" % i, "Nat", rss[i], str)
    E.add("evA", "List (List KLW.WEv)", ll("ev%d" % i for i in range(nA)))
    E.add("evB", "List (List KLW.WEv)", ll("ev%d" % i for i in range(nA, len(segs))))
    E.lst("scatEv", "KLW.WEv", scal[False]["sc_ev"], E.ev, per=150)
    E.lst("wt", "KLW.WOp", list(zip(L.wt, L.wti)), lambda p: "⟨%s, %s⟩" % (ll(map(str, p[0])), ll(map(str, p[1]))), per=80)
    num = str
    E.lst("ports", "Nat", c["ports"], num)
    E.lst("start", "Nat", c["start"], num)
    E.lst("final", "Nat", c["final"], num)
    E.lst("ret", "Nat × Nat", [(s, f) for _, s, f in sorted(c["ret"])], lambda x: "(%d, %d)" % x)
    E.lst("ext", "Nat × Nat", [tuple(x) for x in c["ext"]], lambda x: "(%d, %d)" % x)
    E.lst("scatRows", "List (Nat × GXD.Co)", c["scat"]["table"],
          lambda row: ll("(%d, %s)" % (k, E.cf(co(p, q))) for k, p, q in row), per=60)
    E.add("scat", "GXD.Scat", ".table scatRows")
    E.lst("frames", "List Nat", c["frames"], lambda f: ll(map(str, f)), per=150)
    E.lst("zb", "Nat × Nat", [tuple(x) for x in c["zero_buffers"]], lambda x: "(%d, %d)" % x)
    E.lst("helperOrigin", "Nat", c["helper_origin"], num)
    E.lst("removed", "Nat", c["removed"], num)
    E.lst("levels", "Nat × Nat", sorted((v + int(t), f) for t, f in c["levels"].items()), lambda x: "(%d, %d)" % x)
    E.add("nP", "Nat", str(c["nP"]))
    E.lst("exits", "Nat", sorted(c["exits"]), num)
    # exit-route residuals per slot (as k-combine's make_sched_x / Lean's res_ok): from the start and end label states
    dimS = lambda r: L.dim[states[0][r]]
    resW = [(L.dim[states[-1][n0 - R + q]] if (n0 - R + q) in set(c["exits"]) else h) - dimS(n0 - R + q) for q in range(R)]
    need(resW == c["residuals"], "residuals from the label states differ from the gcert/1w key")
    E.lst("resW", "Nat", resW, num)
    E.lst("unitsEnd", "Nat × GXD.Co", [(r, co(p, q)) for r, p, q in c["units_end"]], lambda x: "(%d, %s)" % (x[0], E.cf(x[1])))
    d = bits(n)
    dF = bits(len(c["frames"]))
    ts = 2
    E.trie("tab", [code(h, f) for f in c["frames"]], dF, low=min(8, dF - ts - 1))
    for i, st in enumerate(states):
        E.trie("st%d" % i, st, d)
    for o, pre in ((False, "rf"), (True, "rb")):
        for i, sn in enumerate(scal[o]["snaps"]):
            rows = [sorted((k, int(x * D)) for k, x in sn[r].items()) for r in range(n)]
            for r in rows:
                for k, z in r:
                    pass
            E.rtrie("%s%d" % (pre, i), rows, d)
    E.add("par", "GLab.Par", "⟨%d, %d, %d, %d, tab⟩" % (h, n, d, h * L.most))
    E.add("spw", "KLW.SPw", "⟨%d, %d, %d, %d, %d⟩" % (d, n, N, B, D))
    # coefficient constants first
    cdefs = ["noncomputable def %s : GXD.Co := ⟨%s, %d, %d⟩" % (nm, "true" if k[0] else "false", k[1], k[2])
             for k, nm in sorted(E.co.items(), key=lambda kv: int(kv[1][1:]))]
    defs = cdefs + E.defs
    mods, cur, size = [], [], 0
    for x in defs:
        if cur and size + len(x) > a.mod_bytes:
            mods.append(cur)
            cur, size = [], 0
        cur.append(x)
        size += len(x)
    mods.append(cur)
    base = os.path.join(a.out, "Work/GCert/Data/Gen")
    os.makedirs(os.path.join(base, NAME), exist_ok=True)
    os.makedirs(os.path.join(a.out, "logs"), exist_ok=True)
    pkg = "Work.GCert.Data.Gen"
    sizes = []
    for i, m in enumerate(mods):
        imp = "import Work.KLowW.Data\nimport Work.KLowW.ScalK\nimport Work.GCert.Labels.Def\n" + \
              ("import %s.%s.D%d\n" % (pkg, NAME, i - 1) if i else "")
        txt = imp + HEAD % dict(tag, what="data segment %d of %d" % (i, len(mods))) + \
            "namespace GXD.%s\nopen GXD\n" % NAME + "\n".join(m) + "\nend GXD.%s\n" % NAME
        open(os.path.join(base, NAME, "D%d.lean" % i), "w").write(txt)
        sizes.append(len(txt.encode()))
    top = "%s.%s.D%d" % (pkg, NAME, len(mods) - 1)
    open(os.path.join(base, NAME + ".lean"), "w").write(
        "import %s\n" % top + HEAD % dict(tag, what="all data (import of the last data segment)") +
        "namespace GXD.%s\nend GXD.%s\n" % (NAME, NAME))
    checks = []

    def chk(mod, imports, doc, thms):
        txt = "".join("import %s\n" % x for x in imports) + HEAD % dict(tag, what=doc) + \
            "namespace GXD.%s\nopen GXD\n" % NAME + "\n".join("theorem %s : %s := by decide +kernel" % t for t in thms) + \
            "\nend GXD.%s\n" % NAME + "".join("#print axioms GXD.%s.%s\n" % (NAME, t[0]) for t in thms)
        open(os.path.join(base, "%s%s.lean" % (NAME, mod)), "w").write(txt)
        checks.append(("%s.%s%s" % (pkg, NAME, mod), [t[0] for t in thms]))
    dat = "%s.%s" % (pkg, NAME)
    for j in range(1 << ts):
        path = bin(j)[2:].zfill(ts)[::-1]
        chk("Tab%d" % j, ["Work.GCert.Labels.Def", dat], "E0: subtrie %d of %d of the frame table" % (j, 1 << ts),
            [("tab%d_ok" % j, "GLab.tabK 9 tab_n%s = true" % path)])
    chk("WT", ["Work.KLowW.LabK", dat], "the warp table: G and G⁻¹ mutually inverse on the unit vectors",
        [("wt_ok", "KLW.wtK 9 wt = true")])
    for i in range(len(segs)):
        chk("Seg%d" % i, ["Work.KLowW.LabK", dat], "label replay of event segment %d (%s, %d events)" % (i, segph[i], len(segs[i])),
            [("seg%d_ok" % i, "KLW.segW par wt ev%d st%d st%d rs%d = true" % (i, i, i + 1, i))])
    for o, pre, ob in ((False, "f", "false"), (True, "b", "true")):
        j = 0
        for i in range(len(segs)):
            if i == nA:
                chk("Sc%sS" % pre.upper(), ["Work.KLowW.ScalK", dat], "scalar replay of the scatter (%s)" % ("dual" if o else "forward"),
                    [("sc%sS_ok" % pre, "KLW.scalW spw %s scatEv r%s%d r%s%d = true" % (ob, pre, j, pre, j + 1))])
                j += 1
            chk("Sc%s%d" % (pre.upper(), i), ["Work.KLowW.ScalK", dat],
                "scalar replay of event segment %d (%s, %s)" % (i, segph[i], "dual" if o else "forward"),
                [("sc%s%d_ok" % (pre, i), "KLW.scalW spw %s ev%d r%s%d r%s%d = true" % (ob, i, pre, j, pre, j + 1))])
            j += 1
    last = len(segs)
    chk("Res", ["Work.Gauge.Labels", "Work.GCert.Labels.Canon", dat],
        "exit-route residuals of the 408 slots from the start and end label states (dim final - dim start for an exit "
        "slot, 9 - dim start otherwise), as T74E's res_ok",
        [("res_ok", "∀ q : Fin %d, Bool.rec (motive := fun _ => Nat) (%d - GLab.cdim (tab.get (st0.get (%d + q.val)))) "
          "(GLab.cdim (tab.get (st%d.get (%d + q.val))) - GLab.cdim (tab.get (st0.get (%d + q.val)))) "
          "(GLab.memK (%d + q.val) exits) = resW.getD q.val 0" % (R, h, n0 - R, last, n0 - R, n0 - R, n0 - R))])
    # ---- histogram facts on the explicit rank lists
    allr = sum(rss, [])
    topr = max(allr)
    lines = ["/-- block ranks of the whole program, in order (the lists of the replay segments) -/",
             "noncomputable def rAll : List Nat := " + " ++ (".join("rs%d" % i for i in range(len(segs))) + ")" * (len(segs) - 1),
             "/-- the dimensions of the 8 centre frames (one c-block each) -/",
             "noncomputable def rC : List Nat := %s" % ll(map(str, cen))]
    thms = [("rankBound", "wsumK (SSC.wGt %d) rAll 0 = 0" % topr)]
    thms += [("rank%d_ok" % r, "wsumK (SSC.wEq %d) rAll 0 = %d" % (r, allr.count(r))) for r in range(topr + 1)]
    thms += [("rankC%d_ok" % r, "wsumK (SSC.wEq %d) rC 0 = %d" % (r, cen.count(r))) for r in sorted(set(cen))]
    txt = "import Work.GCert.Data.Count\nimport %s\n" % dat + HEAD % dict(tag, what="kernel counts of the block ranks "
          "(whole invocation %s)" % dict(sorted(hist.items()))) + "namespace GXD.%s\nopen GXD\n" % NAME + "\n".join(lines) + "\n" + \
          "\n".join("theorem %s : %s := by decide +kernel" % t for t in thms) + "\nend GXD.%s\n" % NAME + \
          "".join("#print axioms GXD.%s.%s\n" % (NAME, t[0]) for t in thms[:2])
    open(os.path.join(base, "%sHist.lean" % NAME), "w").write(txt)
    checks.append(("%s.%sHist" % (pkg, NAME), [t[0] for t in thms]))
    # ---- fingerprint
    txt = "import Work.PWGen.Fp\nimport %s\n" % dat + HEAD % dict(tag, what="the word fingerprint (64-bit FNV-1a, "
          "PWG.fpWord), equal to the Python generator's value") + "namespace GXD.%s\nopen GXD\n" % NAME + \
          "theorem fp_ok : PWG.fpWord evA evB wt ret scatRows start final = %d := by decide +kernel\n" % info["fp"] + \
          "end GXD.%s\n#print axioms GXD.%s.fp_ok\n" % (NAME, NAME)
    open(os.path.join(base, "%sFp.lean" % NAME), "w").write(txt)
    checks.append(("%s.%sFp" % (pkg, NAME), ["fp_ok"]))
    # ---- the label half composed: frame table, warp table, the replays of the two phases (p-warpchk's RunW)
    order = [int(bin(j)[2:].zfill(ts)[::-1], 2) for j in range(1 << ts)]      # Tab<j> checks subtrie path rev(j)
    pos = {bin(j)[2:].zfill(ts)[::-1]: j for j in range(1 << ts)}
    def tnode(path, depth):
        if depth == ts:
            return "tab%d_ok" % pos[path]
        return "(GLab.tabK_node _ _ _ %s %s)" % (tnode(path + "0", depth + 1), tnode(path + "1", depth + 1))
    def nest(xs):
        return " ++ (".join(xs) + ")" * (len(xs) - 1)
    def runs(lo, hi):
        if hi - lo == 1:
            return "(KLW.RunW.of_seg tab_ok wt_ok seg%d_ok)" % lo
        return "(KLW.RunW.append (KLW.RunW.of_seg tab_ok wt_ok seg%d_ok) %s)" % (lo, runs(lo + 1, hi))
    imps = ["Work.KLowW.LabSound", "Work.GCert.Labels.Elim", "%s.%sHist" % (pkg, NAME)] + \
        ["%s.%sTab%d" % (pkg, NAME, j) for j in range(1 << ts)] + ["%s.%sWT" % (pkg, NAME)] + \
        ["%s.%sSeg%d" % (pkg, NAME, i) for i in range(len(segs))]
    ns = len(segs)
    body = f"""theorem tab_ok : GLab.tabK 9 tab = true := {tnode("", 0)}

/-- the events before / after the scatter, flattened -/
theorem evA_eq : evA.flatten = {nest(["ev%d" % i for i in range(nA)])} := by
  show List.flatten [{", ".join("ev%d" % i for i in range(nA))}] = _
  simp only [List.flatten_cons, List.flatten_nil, List.append_nil]
theorem evB_eq : evB.flatten = {nest(["ev%d" % i for i in range(nA, ns)])} := by
  show List.flatten [{", ".join("ev%d" % i for i in range(nA, ns))}] = _
  simp only [List.flatten_cons, List.flatten_nil, List.append_nil]

/-- block ranks of the events before / after the scatter -/
noncomputable def rA : List Nat := {nest(["rs%d" % i for i in range(nA)])}
noncomputable def rB : List Nat := {nest(["rs%d" % i for i in range(nA, ns)])}
theorem rAll_eq : rAll = rA ++ rB := by
  unfold rAll rA rB
  simp only [List.append_assoc]

/-- **the label replay before the scatter** (D0, P, C0, A), composed from its kernel-checked segments -/
theorem runA : KLW.RunW par wt st0 evA.flatten st{nA} rA := by
  rw [evA_eq]
  exact {runs(0, nA)}
/-- **the label replay after the scatter** (B with its warps, tail, final climbs, U with the zero buffers) -/
theorem runB : KLW.RunW par wt st{nA} evB.flatten st{ns} rB := by
  rw [evB_eq]
  exact {runs(nA, ns)}"""
    txt = "".join("import %s\n" % x for x in imps) + HEAD % dict(tag, what="the label half composed (frame table, "
          "warp table, the replays before and after the scatter; p-warpchk's RunW)") + \
        "namespace GXD.%s\nopen GXD\n" % NAME + body + "\nend GXD.%s\n" % NAME + \
        "".join("#print axioms GXD.%s.%s\n" % (NAME, x) for x in ("tab_ok", "runA", "runB"))
    open(os.path.join(base, "%sLabW.lean" % NAME), "w").write(txt)
    checks.append(("%s.%sLabW" % (pkg, NAME), ["tab_ok", "runA", "runB", "rAll_eq"]))
    # ---- the scalar half composed (p-warpchk's scalW_sound / esMap_append): forward and dual, A / scatter / B
    def chain(pre, lo, hi, j0):
        """rewrites for the segments lo..hi-1 starting at row trie j0 (the last one closes the goal)"""
        steps = []
        for t, i in enumerate(range(lo, hi)):
            fact = "(KLW.scalW_sound hq hK sc%s%d_ok).2.2.2 x" % (pre, i)
            steps.append(fact if t == hi - lo - 1 else None)
        rws = []
        for t, i in enumerate(range(lo, hi - 1)):
            rws += ["KLW.esMap_append", "(KLW.scalW_sound hq hK sc%s%d_ok).2.2.2 x" % (pre, i)]
        return rws, "(KLW.scalW_sound hq hK sc%s%d_ok).2.2.2 x" % (pre, hi - 1)
    jB0 = nA + 1
    sc_imps = ["Work.KLowW.ScalSound", "%s.%sLabW" % (pkg, NAME)] + \
        ["%s.%sSc%s%d" % (pkg, NAME, P.upper(), i) for P in "fb" for i in range(ns)] + \
        ["%s.%sSc%sS" % (pkg, NAME, P.upper()) for P in "fb"]
    parts = []
    for pre, o, oname in (("f", "false", "forward"), ("b", "true", "dual")):
        ra, la = chain(pre, 0, nA, 0)
        rb_, lb = chain(pre, nA, ns, jB0)
        R0 = "r%s" % pre
        parts.append(f"""/-- the scalar replay before the scatter ({oname}), from the identity rows `{R0}0` to `{R0}{nA}` -/
theorem scal{pre.upper()}A (hq : KLW.SCtx spw emb un) (hK : ∀ a, K a * K (-a) = 1) {{f : ℕ}} (x : Data α ρ f) :
    KLW.esMap {o} K un evA.flatten f (KLW.vecR spw K un emb {R0}0 x) = KLW.vecR spw K un emb {R0}{nA} x := by
  rw [evA_eq{"".join(", " + r for r in ra)}]
  exact {la}
/-- the scatter ({oname}): `{R0}{nA}` to `{R0}{nA + 1}` -/
theorem scal{pre.upper()}S (hq : KLW.SCtx spw emb un) (hK : ∀ a, K a * K (-a) = 1) {{f : ℕ}} (x : Data α ρ f) :
    KLW.esMap {o} K un scatEv f (KLW.vecR spw K un emb {R0}{nA} x) = KLW.vecR spw K un emb {R0}{nA + 1} x :=
  (KLW.scalW_sound hq hK sc{pre}S_ok).2.2.2 x
/-- the scalar replay after the scatter ({oname}), `{R0}{nA + 1}` to the end rows `{R0}{ns + 1}` -/
theorem scal{pre.upper()}B (hq : KLW.SCtx spw emb un) (hK : ∀ a, K a * K (-a) = 1) {{f : ℕ}} (x : Data α ρ f) :
    KLW.esMap {o} K un evB.flatten f (KLW.vecR spw K un emb {R0}{nA + 1} x) = KLW.vecR spw K un emb {R0}{ns + 1} x := by
  rw [evB_eq{"".join(", " + r for r in rb_)}]
  exact {lb}""")
    facts = f"""/-- the start rows are the identity (all {n} registers), forward and dual -/
theorem rf0_id : ∀ r, r < spw.nr → rf0.get r = [(r, (spw.D : Int))] := by decide +kernel
theorem rb0_id : ∀ r, r < spw.nr → rb0.get r = [(r, (spw.D : Int))] := by decide +kernel
/-- the end rows carry no warp word (every key is a plain register) -/
theorem rf{ns + 1}_free : ∀ r, r < 2 ^ spw.d → ∀ e ∈ rf{ns + 1}.get r, e.1 < spw.N := by decide +kernel
theorem rb{ns + 1}_free : ∀ r, r < 2 ^ spw.d → ∀ e ∈ rb{ns + 1}.get r, e.1 < spw.N := by decide +kernel"""
    txt = "".join("import %s\n" % x for x in sc_imps) + HEAD % dict(tag, what="the scalar half composed (free warp "
          "algebra, forward and dual; p-warpchk's scalW_sound)") + \
        "namespace GXD.%s\nopen SSC GXD GLab OAI.PowerSaving OAI.PowerSaving.Binary OAI.PowerSaving.GF OAI.PowerSaving.RAM\n" % NAME + \
        facts + "\n\nsection\nvariable {α ρ : Type} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]\n" + \
        "variable {K : Int → CMat α} {emb : ρ → Nat} {un : Nat → ρ}\n\n" + "\n\n".join(parts) + "\nend\n" + \
        "end GXD.%s\n" % NAME + "".join("#print axioms GXD.%s.%s\n" % (NAME, x) for x in
                                        ("rf0_id", "scalFA", "scalFS", "scalFB", "scalBA", "scalBS", "scalBB"))
    open(os.path.join(base, "%sScW.lean" % NAME), "w").write(txt)
    checks.append(("%s.%sScW" % (pkg, NAME), ["rf0_id", "rb0_id", "scalFA", "scalFS", "scalFB", "scalBA", "scalBS", "scalBB"]))
    info.update(data_modules=["%s.%s.D%d" % (pkg, NAME, i) for i in range(len(mods))], data_bytes=sizes,
                checks=checks, coefficients=len(E.co))
    json.dump(info, open(os.path.join(a.out, "logs", "wgen-%s.json" % NAME), "w"), indent=1)
    print("wrote %d data modules (%d bytes), %d check modules" % (len(mods), sum(sizes), len(checks)))


if __name__ == "__main__":
    try:
        main()
    except Reject as e:
        print("REJECTED by the mirror (nothing written):", e)
        sys.exit(2)
