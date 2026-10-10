#!/usr/bin/python3 -I
"""build.py: step 1 of the rebuild.  Builds the paired-cube helper word of pull request #168 of the outside
repository with the generator of this directory, and compares it with the graph hash and the profile that the
outside repository records (its files are read as data only).
Parts ADAPTED from scripts of the outside repository (github.com/CrocSwap/integer-mult-bounds, Apache-2.0) at
commit 4a3c769, read as text and never run: target_hist() and the field names of profile() follow the last
part of select() of scripts/paired_cube/gauges.py, and the calls that make the graph and its hash follow
regenerate() of scripts/paired_cube_producer.py, so that the result can be compared with the outside records
field by field.  Notice of both files: "Copyright 2026 icekylinx. Apache-2.0."  NOTICE, section 6.

usage: python3 -I -B build.py pr168        (rebuild.py runs it; see rbpaths.py for the two directories)
writes WORK/cache/pr168.pkl (graph, frames witness, word) for aligned.py.
"""
import sys, os, json, hashlib, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter
from f2 import rref, perp, inside
import graphgen, modules, compile as comp
import rbpaths

PR168 = rbpaths.inputs() + "pr168/"
WORK = rbpaths.work()


def graph_hash(g):
    b = {k: g[k] for k in ("inputs", "labels", "args", "signs", "roots", "centers")}
    return hashlib.sha256(json.dumps(b, separators=(",", ":")).encode()).hexdigest()


def target_hist(g, word, roots):
    h, v, inputs = g["h"], g["v"], g["inputs"]
    Y, cur = Counter(), [()] * v
    for z in reversed(word["selected"]):
        U = perp(z["annihilator"], h)
        for t in z["targets"]:
            assert inside(cur[t], U)
            Y[len(U) - len(cur[t])] += 1
            cur[t] = U
    for r in roots:
        if r["kind"] == "center":
            continue
        U = perp(rref(inputs[t] for t in r["targets"]), h)
        for t in r["targets"]:
            assert inside(cur[t], U)
            Y[len(U) - len(cur[t])] += 1
            cur[t] = U
    for t in range(v):
        Y[h - 1 - len(cur[t])] += 1
    return Y


def profile(g, prof, word, wit):
    h, v, R = g["h"], g["v"], prof["R"]
    Y = target_hist(g, word, wit["roots"])
    C = Counter({r: 3 * c for r, c in word["H"].items() if r})
    C.update({r: 3 * c for r, c in prof["source"].items() if r})
    C.update({r: 3 * c for r, c in Y.items() if r})
    C.update({3 * r: c for r, c in word["gauges"].items()})
    C[2] += 2 * v
    m, W = 3 * h, 2 * v + R
    mass = sum(r * c for r, c in C.items())
    assert W * m - mass == 2 * v - 3 * prof["loss"], (W * m - mass, 2 * v - 3 * prof["loss"])
    touched = set(word["sources"].values())
    for i in word["phase1"]:
        touched.update(word["ops"][i][:2])
    return dict(h=h, v=v, R=R, q=prof["q"], c=prof["c"], matched=prof["matched"], loss=prof["loss"], m=m,
                W_per_vertex=W, rank_per_vertex=mass, deficit_per_vertex=W * m - mass,
                selected_roles=len(word["selected"]),
                selected_rank_histogram={str(k): n for k, n in word["gauges"].items()},
                remaining_internal_histogram=[word["H"][r] for r in range(h + 1)],
                target_data_histogram={str(k): n for k, n in Y.items()},
                source_data_histogram={str(k): n for k, n in prof["source"].items()},
                child_histogram={str(k): n for k, n in C.items()},
                phase1_operations=len(word["phase1"]), total_M_operations=len(word["ops"]),
                untouched_non_source_roles=R - len(touched))


def compare(mine, theirs):
    bad = []
    for k, w in theirs.items():
        if k in ("copied_centers_already", "matching_frames"):
            continue
        a = mine.get(k)
        if isinstance(w, dict):
            a = {kk: n for kk, n in (a or {}).items() if n}
            w = {kk: n for kk, n in w.items() if n}
        if a != w:
            bad.append(k)
    return bad


def build(name):
    if name != "pr168":
        sys.exit("build.py builds the word pr168 only")
    t0 = time.time()
    src = PR168 + "references/paired-cube/sources/"
    G = graphgen.Gen(11)
    G.finish(modules.load_triple(src + "tmod_TE_TD_TB3_1_1_4_full_6.0617964e-4.json", 11),
             modules.load_zero_based(src + "pmod_J0_full_6.0666810e-4.json", 45, modules.pair_want(10)),
             modules.load_zero_based(src + "qmod_climb3u_best.json", 9, [511 ^ (1 << i) for i in range(9)]),
             local=json.load(open(src + "local_L1.json")))
    G.fuse("00111100")
    ref = PR168 + "references/paired-cube/selected-module/"
    arcs = json.load(open(ref + "matching-arcs.json"))
    pin = json.load(open(ref + "SOURCE.json"))["graph_sha256"]
    expected = json.load(open(PR168 + "certificates/paired-cube-complex-input.json"))
    g = G.graph()
    g["kindof"] = G.kindof
    g["counts"] = G.counts
    print("[%s] graph: p=%d h=%d v=%d nodes=%d roots=%d counts=%s  (%.1fs)"
          % (name, g["p"], g["h"], g["v"], len(g["args"]), len(g["roots"]), G.counts, time.time() - t0), flush=True)
    hh = graph_hash(g)
    print("    graph sha256 %s  pinned %s  EQUAL=%s" % (hh[:16], pin[:16], hh == pin), flush=True)
    prof, wit = comp.closure(g, arcs)
    print("    closure: R=%d c=%d q=%d matched=%d loss=%d  (%.1fs)"
          % (prof["R"], prof["c"], prof["q"], prof["matched"], prof["loss"], time.time() - t0), flush=True)
    word = comp.select(g, prof, wit)
    mine = profile(g, prof, word, wit)
    print("    word: ops=%d phase1=%d gauges=%s W=%d rank=%d deficit=%d  (%.1fs)"
          % (len(word["ops"]), len(word["phase1"]), dict(word["gauges"]), mine["W_per_vertex"],
             mine["rank_per_vertex"], mine["deficit_per_vertex"], time.time() - t0), flush=True)
    bad = compare(mine, expected)
    print("    published complex-input profile: %d fields compared, MISMATCHES=%s" % (len(expected), bad), flush=True)
    pickle.dump(dict(g=g, prof=prof, wit=wit, word=word, profile=mine), open(WORK + "cache/%s.pkl" % name, "wb"))
    return g, prof, wit, word, mine


if __name__ == "__main__":
    for nm in sys.argv[1:]:
        build(nm)
