"""stage.py: the literal one-invocation stage word of the paired-cube helper circuit, every register with its
full chain of frames.

Registers: X_t = t, Y_t = v + t, physical slot k = 2v + k (recipients share their donor's slot; sinks have none).
Event kinds (frame F = the common frame at which the gate acts; every register it touches climbs to F first):
  read0   y -= K s            old-value read of an ungauged role at frame 0 (the dirt gate, K = response)
  inj     slot += x_t         at the line <q_t>
  op      dest = ca*dest + cb*ctrl      helper gate, at its operation frame
  centre  y_all += scat * slot          retained star total: slot at the star frame, every y at frame 0
  gread   y_T -= K s          deferred old-value read of a GAUGED role at its entrance gauge (targets T)
  ysh     y_t += sg * y_p     terminal sink shears (y reads y)
  yw      y_p += c * ctrl     terminal sink pivot write
  rread   y_T += (+-1/2) slot           side root delivered at the cap of its target set T
  ktr     X_S := M X_S        4x4 Hadamard/2 on four same-parity ports of a cube, IN PLACE on x roles
  kd      y_t' += x_t         at t'^perp (t' the antipode)
  undo    inverse helper gates and injections at the full frame (no block)
The layer (moved frames, reuse pairs, sinks) is data given by the caller.  The word that rebuild.py builds has
no sinks, so the events ysh and yw do not occur in it.
ADAPTED from two scripts of the outside repository (github.com/CrocSwap/integer-mult-bounds, Apache-2.0) at
commit 4a3c769.  They were read as text and never run.  The order of the events, the backward pass of
responses() and many names are those of research/terminal-sinks/sinks_gate.py (its function sinks_record with
the inner functions build and framescan) and of scripts/paired_cube_physical.py.  Changed here: an event is a
record with named fields; the arithmetic is exact, with fractions, where theirs is modulo a prime; the frame
scan returns every block.  What the two files say of themselves:
  sinks_gate.py            "The substitution is the terminal-sink lemma of PR #166 by jamesyc" and "Gate, audit
                           and integration by eumemic with Anthropic Claude assistance.  Apache-2.0."
  paired_cube_physical.py  "Prepared by eumemic with Anthropic Claude assistance. Apache-2.0. The signed word,
                           frozen arcs, gauges, chronology and shared-core accounting are icekylinx's PR144."
NOTICE, section 6, lists all files of this kind.
"""
from collections import Counter
from fractions import Fraction as Fr
from f2 import rref, perp, inside

ZERO = ()


def responses(g, wit, word):
    """backward pass: resp[s] = (cvec over the h star coordinates or None, {target: coefficient})"""
    h, roots, R = g["h"], wit["roots"], word["R"]
    cseed, dseed = {}, {}
    for r, s in zip(roots, word["rootroles"]):
        if r["kind"] == "center":
            cseed[s] = r["coordinate"]
        else:
            dseed[s] = {t: Fr(q) for t, q in zip(r["targets"], r["coefficients"])}
    cvec, dpart = [None] * R, [dict() for _ in range(R)]
    for s, c in cseed.items():
        cvec[s] = [int(k == c) for k in range(h)]
    for s, d in dseed.items():
        dpart[s] = dict(d)
    for (a, b, _), (ca, cb) in zip(reversed(word["ops"]), reversed(word["coef"])):
        if cvec[a] is not None:
            cvec[b] = [cb * u for u in cvec[a]] if cvec[b] is None else [w + cb * u for w, u in zip(cvec[b], cvec[a])]
        if dpart[a]:
            db = dpart[b]
            for t, u in dpart[a].items():
                db[t] = db.get(t, 0) + cb * u
        if ca != 1:
            if cvec[a] is not None:
                cvec[a] = [ca * u for u in cvec[a]]
            dpart[a] = {t: ca * u for t, u in dpart[a].items()}
    return cseed, dseed, cvec, dpart


def build(g, wit, word, layer, use_sinks=True):
    """returns a dict with the keys events, start, final, nreg, slot, live, merge, sink, gauge, gtargets, frames,
    default, moved, pairs, deadline, position, phase1, rest, resp, rootframe, leaf, FULL, h, v (start and final:
    the frame of every register before the first event and at the end; slot: live role -> register).
    layer = dict(frames=[[op, basis]..], pairs=[[donor, recipient, read_before_op or None]..],
    sinks=[[role, pivot]..])."""
    h, v, R, inputs = g["h"], g["v"], word["R"], g["inputs"]
    FULL = rref(1 << i for i in range(h))
    roots, ann, spans = wit["roots"], wit["ann"], wit["spans"]
    ops, coef = word["ops"], word["coef"]
    frames = [perp(ann[x], h) for _, _, x in ops]
    default = list(frames)
    moved = 0
    for i, F in layer.get("frames", []):
        F = tuple(F)
        assert rref(F) == F and F != frames[i], "moved frame not canonical or not moved"
        frames[i] = F
        moved += 1
    for i, (_, _, x) in enumerate(ops):
        assert inside(spans[x], frames[i]), "value span outside operation frame"
    pairs = [(int(a), int(b)) for a, b, _ in layer.get("pairs", [])]
    deadline = {int(b): (None if t is None else int(t)) for _, b, t in layer.get("pairs", [])}
    merge = {b: a for a, b in pairs}
    phase1 = list(word["phase1"])
    pset = set(phase1)
    rest = [i for i in range(len(ops)) if i not in pset]
    position = {i: k for k, i in enumerate(phase1 + rest)}
    sources = word["sources"]
    leaf = {s: x for x, s in sources.items()}
    selected = [z["role"] for z in word["selected"]]
    gauge = {z["role"]: perp(z["annihilator"], h) for z in word["selected"]}
    gtargets = {z["role"]: z["targets"] for z in word["selected"]}
    deferred = set(selected)
    cseed, dseed, cvec, dpart = responses(g, wit, word)
    rootframe, rootorder = {}, []
    for j, (r, s) in enumerate(zip(roots, word["rootroles"])):
        if r["kind"] == "center":
            rootframe[s] = spans[r["node"]]
        else:
            rootframe[s] = wit["rframe"][j]
            rootorder.append((s, r["targets"]))
    sink = {}
    for s, pv in (layer.get("sinks", []) if use_sinks else []):
        js = [j for j, r in enumerate(word["rootroles"]) if r == s]
        assert len(js) == 1 and roots[js[0]]["kind"] == "side"
        wr = [i for i, (a, b, x) in enumerate(ops) if a == s]
        assert wr and not any(b == s for a, b, x in ops) and all(i not in pset for i in wr)
        assert s not in merge and s not in merge.values() and s not in deferred and s not in leaf and s not in cseed
        assert cvec[s] is None and dpart[s] == dseed[s] and all(coef[i][0] == 1 for i in wr)
        cs = set(roots[js[0]]["coefficients"])
        assert len(cs) == 1 and pv in roots[js[0]]["targets"]
        sink[s] = dict(T=list(roots[js[0]]["targets"]), p=pv, writes=wr, U=rootframe[s], c=Fr(cs.pop()))
    live = [s for s in range(R) if s not in merge and s not in sink]
    slot = {s: 2 * v + i for i, s in enumerate(live)}

    def A(s):
        return slot[merge.get(s, s)]
    ev = []
    for s in range(R):
        if s not in deferred and s not in sink:
            ev.append(dict(k="read0", slot=A(s), role=s, F=ZERO))
    for s in sorted(leaf):
        t = leaf[s] - 1
        ev.append(dict(k="inj", a=A(s), b=t, c=1, F=(inputs[t],), role=s))
    fwd = []

    def op_event(i):
        a, b, x = ops[i]
        fwd.append(i)
        return dict(k="op", a=A(a), b=A(b), ca=coef[i][0], cb=coef[i][1], F=frames[i], op=i, node=x, roles=(a, b))
    for i in phase1:
        ev.append(op_event(i))
    for s, c in cseed.items():
        ev.append(dict(k="centre", slot=A(s), role=s, coord=c, F=rootframe[s], G=ZERO))
    for s, z in sink.items():
        ev += [dict(k="ysh", a=v + t, b=v + z["p"], c=-1, F=ZERO, sink=s) for t in z["T"] if t != z["p"]]
    read_op = {s: (deadline[s] if deadline.get(s) is not None else (rest[0] if rest else None)) for s in selected}
    reads_at = {}
    for s in reversed(selected):
        reads_at.setdefault(read_op[s], []).append(s)
    sink_of_write = {i: s for s, z in sink.items() for i in z["writes"]}
    last_write = {z["writes"][-1]: s for s, z in sink.items()}
    for i in rest:
        for s in reads_at.get(i, ()):
            ev.append(dict(k="gread", slot=A(s), role=s, F=gauge[s], T=gtargets[s]))
        if i in sink_of_write:
            s = sink_of_write[i]
            ev.append(dict(k="yw", a=v + sink[s]["p"], b=A(ops[i][1]), c=coef[i][1] * sink[s]["c"], F=frames[i],
                           op=i, sink=s))
        else:
            ev.append(op_event(i))
        if i in last_write:
            z = sink[last_write[i]]
            ev += [dict(k="ysh", a=v + t, b=v + z["p"], c=1, F=z["U"], sink=last_write[i]) for t in z["T"]
                   if t != z["p"]]
    if not rest:
        for s in reads_at.get(None, ()):
            ev.append(dict(k="gread", slot=A(s), role=s, F=gauge[s], T=gtargets[s]))
    for s, tg in rootorder:
        if s not in sink:
            ev.append(dict(k="rread", slot=A(s), role=s, F=rootframe[s], T=tg))
    for b0 in range(0, v, 8):
        qs = [inputs[b0 + k] for k in range(8)]

        def kcoef(i, j):
            w = bin(qs[i] ^ qs[j]).count("1")
            return 1 if w == 6 else -1 if w == 2 else 0
        for parity in (0, 1):
            S = [k for k in range(8) if bin(k).count("1") % 2 == parity]
            Tg = [k ^ 7 for k in S]
            M = [[Fr(kcoef(j, i), 2) for i in S] for j in Tg]
            regs = [b0 + k for k in S]
            ev.append(dict(k="ktr", regs=regs, M=M, F=rref(qs[k] for k in S)))
            for k in range(4):
                ev.append(dict(k="kd", a=v + b0 + Tg[k], b=b0 + S[k], c=1, F=perp((qs[Tg[k]],), h)))
            ev.append(dict(k="ktr", regs=regs, M=[[M[l][k] for l in range(4)] for k in range(4)], F=FULL))
    for n, i in enumerate(reversed(fwd)):
        a, b, x = ops[i]
        ev.append(dict(k="undo", a=A(a), b=A(b), ca=coef[i][0], cb=coef[i][1], F=FULL, op=i))
    for s in sorted(leaf, reverse=True):
        ev.append(dict(k="uninj", a=A(s), b=leaf[s] - 1, c=-1, F=FULL))
    nreg = 2 * v + len(live)
    start, final = [None] * nreg, [None] * nreg
    for t in range(v):
        start[t], final[t] = (inputs[t],), FULL
        start[v + t], final[v + t] = ZERO, perp((inputs[t],), h)
    for s in live:
        start[slot[s]], final[slot[s]] = gauge.get(s, ZERO), FULL
    return dict(events=ev, start=start, final=final, nreg=nreg, slot=slot, live=live, merge=merge, sink=sink,
                gauge=gauge, gtargets=gtargets, frames=frames, default=default, moved=moved, pairs=pairs,
                deadline=deadline, position=position, phase1=phase1, rest=rest, resp=(cseed, dseed, cvec, dpart),
                rootframe=rootframe, leaf=leaf, FULL=FULL, h=h, v=v)


def touched(e, v):
    """registers of an event with the frame each must stand at: list of (reg, frame)"""
    k = e["k"]
    if k == "read0":
        return [(e["slot"], e["F"])]
    if k in ("gread", "rread"):
        return [(e["slot"], e["F"])] + [(v + t, e["F"]) for t in e["T"]]
    if k == "centre":
        return [(e["slot"], e["F"])]          # every y stands at G = 0: nothing to climb
    if k == "ktr":
        return [(r, e["F"]) for r in e["regs"]]
    return [(e["a"], e["F"]), (e["b"], e["F"])]


def framescan(st):
    """climb every register through its events; returns blocks = [(reg, from, to, event index)], histograms"""
    v, ev = st["v"], st["events"]
    cur = list(st["start"])
    blocks = []
    H = {"x": Counter(), "y": Counter(), "s": Counter(), "c": Counter()}

    def cls(r):
        return "x" if r < v else "y" if r < 2 * v else "s"
    ymoved = False
    for n, e in enumerate(ev):
        for r, F in touched(e, v):
            old = cur[r]
            if old == F:
                continue
            assert inside(old, F), ("frame step not nested", e["k"], r, len(old), len(F))
            blocks.append((r, old, F, n))
            H[cls(r)][len(F) - len(old)] += 1
            cur[r] = F
            ymoved = ymoved or v <= r < 2 * v
        if e["k"] == "centre":
            assert not ymoved, "a y role left frame 0 before the star scatter"
            H["c"][abs(len(e["F"]) - len(e["G"]))] += 1
            blocks.append((-1 - e["coord"], e["G"], e["F"], n))
    for r in range(st["nreg"]):
        F = st["final"][r]
        if cur[r] != F:
            assert inside(cur[r], F), ("final step not nested", r)
            blocks.append((r, cur[r], F, len(ev)))
            H[cls(r)][len(F) - len(cur[r])] += 1
    return blocks, H
