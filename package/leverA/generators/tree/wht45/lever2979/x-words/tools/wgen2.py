#!/usr/bin/env python3
"""p-wgen (round wht45/lever, 19 Oct 2026): a gcert/1w word -> p-warpchk's WARP CERTIFICATE `KLW.WCert` and its kernel
facts `KLW.WValid` (Work/KLowW/Inv.lean), after an exact Python mirror of every kernel function involved.

The whole local program is ONE event list (p-warpchk's design): D0 (twisted adds, expanded), P, C0, A[nP:], the CENTRE
STEPS, B (warps, tail, final climbs), U (cleanup with the 91 zero buffers), then the UNIT SCALES.
* centre step of total k (register s_k at frame F_k): `.za c_k F_k`, `.g (.out F_k s_k [(c_k, 1)] [])`, `.dn c_k 0`
  (one block of rank dim F_k: the c-block), `.g (.out 0 c_k [(y_t, coef) for every target t of k] [])`, `.zf c_k`;
  centre copies are the registers 723..730 (after the zero buffers 632..722);
* unit scales: gcert/1w leaves 42 exit slots with unit u != 1 (`units_end`); `.sc f r <u>` at the slot's final frame
  (no climb, free) makes every slot row the identity, so the end rows are p-warpchk's shear rows `shRow`.
Registers: 0..111 x, 112..223 y, 224..631 slots, 632..722 zero buffers, 723..730 centre copies (nc = 99, n = 731).

The mirror classes are tools/wgen.py's (label replay incl. warp transport; free-warp-algebra replay) with `.dn` and
`.sc` added from p-warpchk's LabK.lean / ScalK.lean.  Emits the data, one kernel-check module per WValid field (plus
the 2,000-event segment checks as a fallback), the histogram, residual and fingerprint modules, and `<N>WC`
(`wc : KLW.WCert`, `wvalid : KLW.WValid wc`).

usage: wgen2.py <word.gcert1w.json.gz> <NAME> <out dir> [--seg N] [--mirror-only] [--expect JSON] [--mod-bytes B]"""
import sys, os, json, gzip, hashlib, argparse
from fractions import Fraction as Q
from math import lcm
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wgen as W
from wgen import need, Reject, co, ll, lint, bits, code, gate_adds


class Lab2(W.Lab):
    def event(s, e, S, rs):
        k = e[0]
        if k == "sc":
            s.gate(("out", e[1], e[2], [], []), S, rs)
        elif k == "cut":
            _, r, f, tg = e
            need(r < s.n and S[r] == f and all(t < s.n and S[t] == 0 for t, _ in tg), ("cut", r, f))
            rs.append(s.dim[f])
        elif k == "dn":
            _, r, f = e
            need(r < s.n, "dn register")
            cur = S[r]
            need(s.dim[f] < s.dim[cur] and all(W.in_span(s.fr[cur], u) for u in s.fr[f]), ("dn not inside", e))
            rs.append(s.dim[cur] - s.dim[f])
            S[r] = f
        else:
            W.Lab.event(s, e, S, rs)


class Scal2(W.Scal):
    def event(s, e, R, o):
        k = e[0]
        if k == "dn":
            return
        if k == "cut":
            for t, cf in e[3]:
                if not o:
                    s.add(R, t, e[1], cf)
                else:
                    s.add(R, e[1], t, (not cf[0], cf[1], cf[2]))
            return
        if k == "sc":
            _, f, r, cf = e
            ng, num, den = cf
            need(r < s.nr and num != 0 and den != 0, "sc")
            R[r] = s.scaled(R[r], ng, num, den) if not o else s.scaled(R[r], ng, den, num)
            return
        W.Scal.event(s, e, R, o)


class Emit2(W.Emit):
    def ev(s, e):
        if e[0] == "sc":
            return ".sc %d %d %s" % (e[1], e[2], s.cf(e[3]))
        if e[0] == "dn":
            return ".dn %d %d" % (e[1], e[2])
        if e[0] == "cut":
            return ".cut %d %d %s" % (e[1], e[2], ll("(%d, %s)" % (t, s.cf(k)) for t, k in e[3]))
        return W.Emit.ev(s, e)


def fp_nats(evs, wt, wti, c):
    """PWG.fpWord's stream: events, then wt (G, Gi), ret, scat rows, start, final"""
    out = []
    def co_(k):
        out.extend([1 if k[0] else 0, k[1], k[2]])
    def lst(xs, f):
        out.append(len(xs))
        for x in xs:
            f(x)
    def gate(g):
        if g[0] == "out":
            out.extend([0, g[1], g[2]]); lst(g[3], lambda p: (out.append(p[0]), co_(p[1]))); lst(g[4], out.append)
        else:
            out.extend([1, g[1], g[2]]); lst(g[3], lambda p: (out.append(p[0]), co_(p[1])))
    def ev(e):
        k = e[0]
        if k == "g":
            out.append(0); gate(e[1])
        elif k == "sc":
            out.extend([5, e[1], e[2]]); co_(e[3])
        elif k == "w":
            out.extend([2, e[1], W.enc(e[2]), e[3], e[4]])
        elif k == "za":
            out.extend([3, e[1], e[2]])
        elif k == "zf":
            out.extend([4, e[1]])
        elif k == "dn":
            out.extend([6, e[1], e[2]])
        else:
            out.extend([7, e[1], e[2]]); lst(e[3], lambda p: (out.append(p[0]), co_(p[1])))
    lst(evs, ev)
    lst(list(zip(wt, wti)), lambda p: (lst(p[0], out.append), lst(p[1], out.append)))
    lst(sorted(c["ret"]), lambda r: out.extend([r[1], r[2]]))
    lst(c["scat"]["table"], lambda row: lst(row, lambda x: (out.append(x[0]), co_(co(x[1], x[2])))))
    lst(c["start"], out.append)
    lst(c["final"], out.append)
    return out


def shrow(D, v, o, r):
    """p-warpchk's shRow (Inv.lean): forward y += x, backward x -= y; slots identity"""
    if r < v:
        return [(r, D)] if not o else [(r, D), (r + v, -D)]
    if r < 2 * v:
        return [(r - v, D), (r, D)] if not o else [(r, D)]
    return [(r, D)]


HEAD = W.HEAD


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cert")
    ap.add_argument("name")
    ap.add_argument("out")
    ap.add_argument("--seg", type=int, default=2000)
    ap.add_argument("--mod-bytes", type=int, default=700000)
    ap.add_argument("--mirror-only", action="store_true")
    ap.add_argument("--test-no-dual", action="store_true",
                    help="TEST ONLY: skip the dual end-row assertion and emit no dual end facts and no WC module")
    ap.add_argument("--centre", choices=["copy", "cut"], default="copy",
                    help="copy: p-warpchk's spec (forward copies .za/.dn/.zf in evs, backward copies in evsB); "
                         "cut: one `.cut s f targets` event per centre (an alternative, not in p-warpchk's Data)")
    ap.add_argument("--expect", required=True, help="expected local profile (JSON); x-words: required, no #3945 default")
    a = ap.parse_args()
    raw = open(a.cert, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    c = json.loads(gzip.decompress(raw))
    need(c.get("format") == "gcert/1w" and c["copies"] == 0, "not a gcert/1w file without copies")
    NAME = a.name
    h, v, R, nZ = c["h"], c["v"], c["R"], c["nZ"]
    n0 = 2 * v + R
    nC = len(c["ret"])
    nc = nZ + nC
    n = n0 + nc
    need([z for z, _ in c["zero_buffers"]] == list(range(n0, n0 + nZ)), "zero buffer numbering")
    phases = dict(W.events(c))
    levels = {int(t): f for t, f in c["levels"].items()}
    # ---- centre steps
    ret = sorted(c["ret"])
    need([k for k, _, _ in ret] == list(range(nC)), "retained totals numbered 0..")
    tgts = {k: [] for k in range(nC)}
    for t, row in enumerate(c["scat"]["table"]):
        for k, p, q in row:
            need(levels.get(t, 0) == 0, ("a centre serves a lifted target", t))
            tgts[k].append((t, co(p, q)))
    CS, CSB = [], []
    for k, s, f in (ret if a.centre == "copy" else []):
        # BACKWARD centre step (its DUAL reading is the hub's backward step): erase the copy at frame 0, gather from the
        # targets (dual of y_t += a c is c += -a y_t), climb UP to F (the rank-8 block), subtract into the total (dual of
        # c += -1 s is s += c), erase the copy back at frame 0
        ck = n0 + nZ + k
        CSB += [("za", ck, 0)]
        if tgts[k]:
            CSB += [("g", ("out", 0, ck, [(v + t, cf) for t, cf in tgts[k]], []))]
        CSB += [("g", ("out", f, ck, [], [])), ("g", ("out", f, s, [(ck, (True, 1, 1))], [])), ("za", ck, 0)]
        ck = n0 + nZ + k
        CS += [("za", ck, f), ("g", ("out", f, s, [(ck, (False, 1, 1))], [])), ("dn", ck, 0)]
        if tgts[k]:
            CS += [("g", ("out", 0, ck, [(v + t, cf) for t, cf in tgts[k]], []))]
        CS += [("zf", ck)]
    if a.centre == "cut":
        CS = [("cut", s, f, [(v + t, cf) for t, cf in tgts[k]]) for k, s, f in ret]
        CSB = CS
        nc = nZ                                   # no centre copies
        n = n0 + nc
    final = c["final"]
    US = [("sc", final[r], r, co(p, q)) for r, p, q in c["units_end"]]
    prog = [("D0", phases["D0"]), ("P", phases["P"]), ("C0", phases["C0"]), ("A", phases["A"]), ("CS", CS),
            ("B", phases["B"]), ("U", phases["U"]), ("US", US)]
    segs, segph = [], []
    for nm, evs in prog:
        for i in range(0, len(evs), a.seg):
            segs.append(evs[i:i + a.seg])
            segph.append(nm)
    allev = sum(segs, [])
    need(len(CS) <= a.seg and len(CSB) <= a.seg, "centre steps longer than a segment")
    iCS = segph.index("CS") if CS else None
    segsB = [CSB if i == iCS else sg for i, sg in enumerate(segs)]      # the backward program (centre steps differ)
    allevB = sum(segsB, [])
    # ---- label mirror
    L = Lab2(c, n)
    S = list(c["start"]) + [0] * nc
    states, rss = [list(S)], []
    for sg in segs:
        rs = []
        for e in sg:
            L.event(e, S, rs)
        rss.append(rs)
        states.append(list(S))
    need(states[-1][:n0] == final, "the label replay does not end at the final frames")
    SB = list(c["start"]) + [0] * nc
    statesB, rssB = [list(SB)], []
    for sg in segsB:
        rs = []
        for e in sg:
            L.event(e, SB, rs)
        rssB.append(rs)
        statesB.append(list(SB))
    need(statesB == states and rssB == rss, "the backward program does not have the forward's states and ranks at the cuts")
    allr = sum(rss, [])
    hist = Counter(allr)
    need({str(r): x for r, x in sorted(hist.items())} == c["hub_local_profile"], ("histogram", dict(hist)))
    need(dict(sorted(hist.items())) == {int(k): x for k, x in json.loads(a.expect).items()}, ("profile", a.expect))
    # ---- scalar mirror, exact, both orientations
    N, B = n, 256
    snaps, Ds = {}, []
    for o in (False, True):
        Sc = Scal2(n, N, B)
        Rw = [{r: Q(1)} for r in range(n)]
        sn = [[dict(x) for x in Rw]]
        for sg in (segsB if o else segs):
            for e in sg:
                Sc.event(e, Rw, o)
            sn.append([dict(x) for x in Rw])
        snaps[o] = sn
        Ds.append(Sc.dens)
    D = lcm(*Ds)
    for o in (False, True):
        end = snaps[o][-1]
        for r in range(n):
            row = sorted((k, int(x * D)) for k, x in end[r].items())
            need(all(x * D == int(x * D) for x in end[r].values()), "end row not integral")
            if o and a.test_no_dual:
                continue
            if r < n0:
                need(row == shrow(D, v, o, r), ("end row is not the shear row", o, r, row[:4]))
            else:
                need(row == [], ("scratch row not empty at the end", o, r))
            need(all(k < N for k, _ in row), "worded end row")
    gl = [r for r, _ in c["ext"]]
    el = sorted(c["exits"])
    dimS = lambda r: L.dim[states[0][r]]
    resW = [(L.dim[states[-1][n0 - R + q]] if (n0 - R + q) in set(el) else h) - dimS(n0 - R + q) for q in range(R)]
    need(resW == c["residuals"], "residuals from the label states differ from the gcert/1w key")
    info = dict(name=NAME, cert=a.cert, cert_sha256=sha, provider=c["provider"], n=n, nc=nc, events=len(allev),
                segments=len(segs), seg_phases=segph, seg_lengths=[len(x) for x in segs],
                ranks=dict(sorted(hist.items())), D=D, B=B, N=N, K=h * L.most,
                max_row_entries=[max(sum(len(x) for x in s_) for s_ in snaps[o]) for o in (False, True)])
    ns_ = fp_nats(allev + allevB, L.wt, L.wti, c)
    info["fp"] = W.fnv(ns_)
    info["fp_stream_length"] = len(ns_)
    print(json.dumps(info, indent=1))
    if a.mirror_only:
        return
    # ---- Lean data
    E = Emit2()
    tag = dict(name=NAME, prov=c["provider"], src=os.path.basename(a.cert), sha=sha)
    num = str
    nseg = len(segs)
    for i, sg in enumerate(segs):
        E.lst("ev%d" % i, "KLW.WEv", sg, E.ev, per=150)
        E.lst("rs%d" % i, "Nat", rss[i], num)
    E.add("evs", "List KLW.WEv", " ++ (".join("ev%d" % i for i in range(nseg)) + ")" * (nseg - 1))
    segB_names = ["ev%d" % i if i != iCS else "evCSB" for i in range(nseg)]
    if iCS is not None:
        E.lst("evCSB", "KLW.WEv", CSB, E.ev, per=150)
    E.add("evsB", "List KLW.WEv", " ++ (".join(segB_names) + ")" * (nseg - 1))
    E.add("rAll", "List Nat", " ++ (".join("rs%d" % i for i in range(nseg)) + ")" * (nseg - 1))
    E.lst("wt", "KLW.WOp", list(zip(L.wt, L.wti)), lambda p: "⟨%s, %s⟩" % (ll(map(str, p[0])), ll(map(str, p[1]))), per=80)
    E.lst("ports", "Nat", c["ports"], num)
    E.lst("start", "Nat", c["start"], num)
    E.lst("final", "Nat", final, num)
    E.lst("ret", "Nat × Nat", [(s, f) for _, s, f in ret], lambda x: "(%d, %d)" % x)
    E.lst("ext", "Nat × Nat", [tuple(x) for x in c["ext"]], lambda x: "(%d, %d)" % x)
    E.lst("scatRows", "List (Nat × GXD.Co)", c["scat"]["table"],
          lambda row: ll("(%d, %s)" % (k, E.cf(co(p, q))) for k, p, q in row), per=60)
    E.lst("frames", "List Nat", c["frames"], lambda f: ll(map(str, f)), per=150)
    E.lst("zb", "Nat × Nat", [tuple(x) for x in c["zero_buffers"]], lambda x: "(%d, %d)" % x)
    if a.centre == "copy":
        E.lst("centreCopies", "Nat × Nat", [(n0 + nZ + k, s) for k, s, _ in ret], lambda x: "(%d, %d)" % x)
    E.lst("helperOrigin", "Nat", c["helper_origin"], num)
    E.lst("removed", "Nat", c["removed"], num)
    E.lst("levels", "Nat × Nat", sorted((v + t, f) for t, f in levels.items()), lambda x: "(%d, %d)" % x)
    E.add("nP", "Nat", str(c["nP"]))
    E.lst("gl", "Nat", gl, num)
    E.lst("exits", "Nat", el, num)
    E.lst("resW", "Nat", resW, num)
    E.lst("unitsEnd", "Nat × GXD.Co", [(r, co(p, q)) for r, p, q in c["units_end"]], lambda x: "(%d, %s)" % (x[0], E.cf(x[1])))
    d = bits(n)
    dF = bits(len(c["frames"]))
    ts = 2
    E.trie("tab", [code(h, f) for f in c["frames"]], dF, low=min(8, dF - ts - 1))
    for i, st in enumerate(states):
        E.trie("st%d" % i, st, d)
    for o, pre in ((False, "rf"), (True, "rb")):
        for i, sn in enumerate(snaps[o]):
            if o and i == 0:
                continue                                  # rb0 = rf0 (the identity rows)
            rows = [sorted((k, int(x * D)) for k, x in sn[r].items()) for r in range(n)]
            E.rtrie("%s%d" % (pre, i), rows, d)
    E.add("rb0", "KLW.RT", "rf0")
    E.add("par", "GLab.Par", "⟨%d, %d, %d, %d, tab⟩" % (h, n, d, h * L.most))
    E.add("pe", "GLab.Par", "⟨%d, %d, %d, %d, tab⟩" % (h, n0, d, h * L.most))
    E.add("spw", "KLW.SPw", "⟨%d, %d, %d, %d, %d⟩" % (d, n, N, B, D))
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
    open(os.path.join(base, NAME + ".lean"), "w").write(
        "import %s.%s.D%d\n" % (pkg, NAME, len(mods) - 1) + HEAD % dict(tag, what="all data") +
        "namespace GXD.%s\nend GXD.%s\n" % (NAME, NAME))
    dat = "%s.%s" % (pkg, NAME)
    checks = []

    def chk(mod, imports, doc, thms):
        txt = "".join("import %s\n" % x for x in imports) + HEAD % dict(tag, what=doc) + \
            "namespace GXD.%s\nopen GXD\n" % NAME + "\n".join("theorem %s : %s := by decide +kernel" % t for t in thms) + \
            "\nend GXD.%s\n" % NAME + "".join("#print axioms GXD.%s.%s\n" % (NAME, t[0]) for t in thms)
        open(os.path.join(base, "%s%s.lean" % (NAME, mod)), "w").write(txt)
        checks.append(("%s.%s%s" % (pkg, NAME, mod), [t[0] for t in thms]))
    last = nseg
    for j in range(1 << ts):
        path = bin(j)[2:].zfill(ts)[::-1]
        chk("Tab%d" % j, ["Work.GCert.Labels.Def", dat], "E0: subtrie %d of %d of the frame table" % (j, 1 << ts),
            [("tab%d_ok" % j, "GLab.tabK 9 tab_n%s = true" % path)])
    chk("WT", ["Work.KLowW.LabK", dat], "the warp table", [("wt_ok", "KLW.wtK 9 wt = true")])
    # WValid fields: one kernel evaluation each over the WHOLE event list
    chk("Lab", ["Work.KLowW.LabK", dat], "THE label replay of the whole event list (%d events)" % len(allev),
        [("lab_ok", "KLW.segW par wt evs st0 st%d rAll = true" % last)])
    chk("LabB", ["Work.KLowW.LabK", dat], "THE label replay of the whole BACKWARD event list (%d events)" % len(allevB),
        [("labB_ok", "KLW.segW par wt evsB st0 st%d rAll = true" % last)])
    chk("ScF", ["Work.KLowW.ScalK", dat], "THE free-warp-algebra replay of the whole event list, forward",
        [("scF_ok", "KLW.scalW spw false evs rf0 rf%d = true" % last)])
    if not a.test_no_dual:
        chk("ScB", ["Work.KLowW.ScalK", dat], "THE free-warp-algebra replay of the whole BACKWARD event list, dual",
            [("scB_ok", "KLW.scalW spw true evsB rb0 rb%d = true" % last)])
    chk("Ends", ["Work.KLow.LabelsX", dat], "the exit-aware end check on the x, y and slot registers",
        [("ends_ok", "GLab.endKX pe %d ports [] gl exits st0 st0 st%d = true" % (v, last))])
    chk("Rows", ["Work.KLowW.Inv", dat], "start rows identity; end rows well formed, shear rows, warp-free",
        [("id0_ok", "KLW.idK spw rf0 = true"), ("wfF_ok", "KLW.RT.wfK rf%d spw.d = true" % last),
         ("shF_ok", "KLW.shKn spw %d false rf%d %d = true" % (v, last, n0)),
         ("frF_ok", "KLW.freeK spw rf%d = true" % last)] +
        ([] if a.test_no_dual else [("wfB_ok", "KLW.RT.wfK rb%d spw.d = true" % last),
                                    ("shB_ok", "KLW.shKn spw %d true rb%d %d = true" % (v, last, n0)),
                                    ("frB_ok", "KLW.freeK spw rb%d = true" % last)]))
    # fallback: the same replays in 2,000-event segments (composable with RunW.append / esMap_append)
    for i in range(nseg):
        chk("Seg%d" % i, ["Work.KLowW.LabK", dat], "label replay of segment %d (%s, %d events)" % (i, segph[i], len(segs[i])),
            [("seg%d_ok" % i, "KLW.segW par wt ev%d st%d st%d rs%d = true" % (i, i, i + 1, i))])
    for pre, ob in (("f", "false"), ("b", "true")):
        for i in range(nseg):
            chk("Sc%s%d" % (pre.upper(), i), ["Work.KLowW.ScalK", dat], "scalar replay of segment %d (%s)" % (i, segph[i]),
                [("sc%s%d_ok" % (pre, i), "KLW.scalW spw %s %s r%s%d r%s%d = true"
                  % (ob, segB_names[i] if pre == "b" else "ev%d" % i, pre, i, pre, i + 1))])
    chk("Res", ["Work.Gauge.Labels", "Work.GCert.Labels.Canon", dat], "exit-route residuals of the 408 slots",
        [("res_ok", "∀ q : Fin %d, Bool.rec (motive := fun _ => Nat) (%d - GLab.cdim (tab.get (st0.get (%d + q.val)))) "
          "(GLab.cdim (tab.get (st%d.get (%d + q.val))) - GLab.cdim (tab.get (st0.get (%d + q.val)))) "
          "(GLab.memK (%d + q.val) exits) = resW.getD q.val 0" % (R, h, n0 - R, last, n0 - R, n0 - R, n0 - R))])
    topr = max(allr)
    thms = [("rankBound", "wsumK (SSC.wGt %d) rAll 0 = 0" % topr)]
    thms += [("rank%d_ok" % r, "wsumK (SSC.wEq %d) rAll 0 = %d" % (r, allr.count(r))) for r in range(topr + 1)]
    chk("Hist", ["Work.GCert.Data.Count", dat], "kernel counts of the block ranks (whole invocation %s, the 8 centre "
        "blocks included)" % dict(sorted(hist.items())), thms)
    txt = "import Work.PWGen.Fp\nimport %s\n" % dat + HEAD % dict(tag, what="the word fingerprint (64-bit FNV-1a, "
          "PWG.fpWord), equal to the Python generator's value") + "namespace GXD.%s\nopen GXD\n" % NAME + \
          "theorem fp_ok : PWG.fpWord (evs ++ evsB) wt ret scatRows start final = %d := by decide +kernel\n" % info["fp"] + \
          "end GXD.%s\n#print axioms GXD.%s.fp_ok\n" % (NAME, NAME)
    open(os.path.join(base, "%sFp.lean" % NAME), "w").write(txt)
    checks.append(("%s.%sFp" % (pkg, NAME), ["fp_ok"]))
    # ---- the certificate and its kernel facts
    order = {"0": 0}
    pos = {bin(j)[2:].zfill(ts)[::-1]: j for j in range(1 << ts)}
    def tnode(path, depth):
        if depth == ts:
            return "tab%d_ok" % pos[path]
        return "(GLab.tabK_node _ _ _ %s %s)" % (tnode(path + "0", depth + 1), tnode(path + "1", depth + 1))
    imps = ["Work.KLowW.Inv", "Work.GCert.Labels.Elim"] + ["%s.%s%s" % (pkg, NAME, m) for m in
            ["Tab0", "Tab1", "Tab2", "Tab3", "WT", "Lab", "LabB", "ScF", "ScB", "Ends", "Rows", "Hist", "Res", "Fp"]]
    body = f"""theorem tab_ok : GLab.tabK 9 tab = true := {tnode("", 0)}

/-- **the warp certificate** of the word (p-warpchk's `KLW.WCert`) -/
noncomputable def wc : KLW.WCert :=
  {{ v := {v}, R := {R}, nc := {nc}, p := par, wt := wt, evs := evs, evsB := evsB, S0 := st0, SF := st{last}, rs := rAll,
    ports := ports, gl := gl, el := exits, q := spw, R0 := rf0, RF := rf{last}, RB := rb{last} }}

/-- **its kernel facts** (every field a `decide +kernel` check of its own module) -/
theorem wvalid : KLW.WValid wc where
  hn := rfl
  nr := rfl
  tab := tab_ok
  wtab := wt_ok
  lab := lab_ok
  labB := labB_ok
  ends := ends_ok
  scF := scF_ok
  scB := scB_ok
  id0 := id0_ok
  wfF := wfF_ok
  wfB := wfB_ok
  shF := shF_ok
  shB := shB_ok
  frF := frF_ok
  frB := frB_ok"""
    txt = "".join("import %s\n" % x for x in imps) + HEAD % dict(tag, what="the warp certificate `wc : KLW.WCert` and "
          "its kernel facts `wvalid : KLW.WValid wc`") + "namespace GXD.%s\nopen GXD\n" % NAME + body + \
          "\nend GXD.%s\n#print axioms GXD.%s.wvalid\n" % (NAME, NAME)
    if not a.test_no_dual:
        open(os.path.join(base, "%sWC.lean" % NAME), "w").write(txt)
        checks.append(("%s.%sWC" % (pkg, NAME), ["tab_ok", "wvalid"]))
    info.update(data_modules=["%s.%s.D%d" % (pkg, NAME, i) for i in range(len(mods))], data_bytes=sizes,
                checks=checks, coefficients=len(E.co), last=last)
    json.dump(info, open(os.path.join(a.out, "logs", "wgen-%s.json" % NAME), "w"), indent=1)
    print("wrote %d data modules (%d bytes), %d check modules" % (len(mods), sum(sizes), len(checks)))


if __name__ == "__main__":
    try:
        main()
    except Reject as e:
        print("REJECTED by the mirror (nothing written):", e)
        sys.exit(2)
