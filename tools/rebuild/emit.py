#!/usr/bin/python3 -I
"""emit.py: step 5 of the rebuild.  Writes the circuit as a JSON file of format gcert/0, the format that
tools/gx/gxconv.py converts to gcert/1.  gcert/0 extends carrier-cert/1 (the format that tools/ccheck.py checks):
labels are subspace BASES, slots may start at an entrance gauge, slots may be reused, gates into x and gates
reading y are explicit event kinds.

usage: python3 -I -B emit.py <name> <out.json>        (after phys.py; <name> becomes the field `name` of the file)
"""
import sys, os, json, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from f2 import kind, radical, step_dirs
import stage
import rbpaths

WORK = rbpaths.work()


def q(c):
    c = Fr(c)
    return [c.numerator, c.denominator]


def emit(name, out):
    t0 = time.time()
    d = pickle.load(open(WORK + "cache/%s.pkl" % name, "rb"))
    g, prof, wit, word = d["g"], d["prof"], d["wit"], d["word"]
    sp = pickle.load(open(WORK + "cache/%s.stage.pkl" % name, "rb"))
    layer = sp["layer"]
    st = stage.build(g, wit, word, layer)
    blocks, H = stage.framescan(st)
    h, v = g["h"], g["v"]
    cseed, dseed, cvec, dpart = st["resp"]
    ftab, fid = [(), st["FULL"]], {(): 0, st["FULL"]: 1}

    def F(fr):
        fr = tuple(fr)
        if fr not in fid:
            fid[fr] = len(ftab)
            ftab.append(fr)
        return fid[fr]
    ev = []
    for e in st["events"]:
        k = e["k"]
        if k == "read0":
            s = e["role"]
            cp = [] if cvec[s] is None else [[i, u] for i, u in enumerate(cvec[s]) if u]
            dp = [[t] + q(u) for t, u in sorted(dpart[s].items()) if u]
            ev.append(["read0", e["slot"], cp, dp])
        elif k == "inj":
            ev.append(["inj", e["a"], e["b"], F(e["F"])])
        elif k == "op":
            ev.append(["op", e["a"], e["b"], e["ca"], e["cb"], F(e["F"]), g["kindof"].get(e["node"] - 1, "port")])
        elif k == "centre":
            ev.append(["centre", e["slot"], e["coord"], F(e["F"])])
        elif k == "gread":
            ev.append(["gread", e["slot"], F(e["F"]), [[t] + q(-u) for t, u in sorted(dpart[e["role"]].items()) if u],
                       list(e["T"])])
        elif k == "rread":
            ev.append(["rread", e["slot"], F(e["F"]), [[t] + q(u) for t, u in sorted(dseed[e["role"]].items())]])
        elif k == "ysh":
            ev.append(["ysh", e["a"], e["b"], e["c"], F(e["F"])])
        elif k == "yw":
            ev.append(["yw", e["a"], e["b"], q(e["c"]), F(e["F"])])
        elif k == "ktr":
            ev.append(["ktr", list(e["regs"]), [[q(c) for c in row] for row in e["M"]], F(e["F"])])
        elif k == "kd":
            ev.append(["kd", e["a"], e["b"], F(e["F"])])
        elif k == "undo":
            ev.append(["undo", e["a"], e["b"], e["ca"], e["cb"]])
        elif k == "uninj":
            ev.append(["uninj", e["a"], e["b"]])
        else:
            raise ValueError(k)
    chains, memo = [], {}
    for r, A, B, n in blocks:          # every block with its new directions (legal ones: orthonormal)
        if (A, B) not in memo:
            memo[A, B] = step_dirs(A, B, h)
        kd, dirs = memo[A, B]
        chains.append([r, n, F(A), F(B), kd, dirs])
    start = [F(x) for x in st["start"]]
    final = [F(x) for x in st["final"]]
    hosts = {}
    for s in st["live"]:
        hosts[st["slot"][s]] = [s]
    for b, a in st["merge"].items():
        hosts[st["slot"][a]].append(b)
    gread_at = {e["role"]: n for n, e in enumerate(st["events"]) if e["k"] == "gread"}
    gauge = []
    for z in word["selected"]:
        s = z["role"]
        slot = st["slot"][st["merge"].get(s, s)]
        gauge.append([s, F(st["gauge"][s]), list(z["targets"]), slot, gread_at.get(s, -1), s not in st["merge"]])
    reuse = [[st["slot"][a], a, b, gread_at[b], -1 if st["deadline"].get(b) is None else st["deadline"][b]]
             for a, b in st["pairs"]]
    ret = [[e["coord"], e["slot"], F(e["F"])] for e in st["events"] if e["k"] == "centre"]
    fk = []
    for fr in ftab:
        fk.append([kind(fr, h), len(radical(fr, h))])
    S = sp["summary"]["with sinks" if "with sinks" in sp["summary"] else "no sinks"]
    cert = dict(
        format="gcert/0", extends="carrier-cert/1", name=name, p=g["p"], h=h, v=v, R=len(st["live"]),
        R_word=word["R"], cst=prof["loss"], N=S["N_moves"],
        registers="x_t = t (t < v); y_t = v + t; slot k = 2v + k (k < R); vectors are bit masks, coordinate c = bit c",
        ports=g["inputs"], frames=[list(fr) for fr in ftab], frame_kind=fk, start=start, final=final,
        scat=dict(inside=[1, 3], outside=[-1, 6],
                  rule="event centre k: every y_t += (1/3 if coordinate k in port t else -1/6) * slot"),
        events=ev, hosts=[hosts[2 * v + i] for i in range(len(st["live"]))], gauge=gauge, reuse=reuse, ret=ret,
        sinks=[[s, z["p"], list(z["T"])] for s, z in st["sink"].items()],
        marks=dict(into_x="events of kind ktr (4x4 gate on four x roles)", from_y="events of kind ysh (y reads y)",
                   slot_to_y="gread, rread, yw", x_to_y="kd", x_to_slot="inj",
                   exterior_gauges=sum(1 for z in gauge if z[5])),
        chains=chains, blocks=S["H"], gauge_tails=S["gauge_tails"],
        counts=dict(events=len(ev), frames=len(ftab), moved_frames=st["moved"], pairs=len(st["pairs"]),
                    sinks=len(st["sink"]), gauged_roles=len(gauge), W_three_stage=S["W_per_vertex"],
                    deficit_three_stage=S["deficit_per_vertex"], rank_three_stage=S["rank_per_vertex"],
                    event_kinds={k: sum(1 for e in ev if e[0] == k) for k in sorted(set(e[0] for e in ev))}),
        classification=sp["summary"].get("classification"),
        info=dict(generator="tools/rebuild (graphgen.py, compile.py, aligned.py, rewrite.py, stage.py, emit.py)",
                  outside_data="data files of CrocSwap/integer-mult-bounds, loaded as data; tools/gx/ORIGIN.md lists them",
                  graph_counts=g["counts"], nodes=len(g["args"]), layer_note=layer.get("note", "")))
    json.dump(cert, open(out, "w"), separators=(",", ":"))
    print("[%s] wrote %s: %d bytes, %d events, %d frames, R=%d  (%.1fs)"
          % (name, out, os.path.getsize(out), len(ev), len(ftab), cert["R"], time.time() - t0), flush=True)


if __name__ == "__main__":
    emit(*sys.argv[1:3])
