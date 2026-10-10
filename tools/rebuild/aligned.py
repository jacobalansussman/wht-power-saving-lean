#!/usr/bin/python3 -I
"""aligned.py: step 2 of the rebuild.  The "source-parity local word" of pull requests #184 / #191 / #193 of the
outside repository on the modules of #168: local circuit with all G channels 's', all A channels 'fd', F = 2;
carrier arcs and preferred operation frames transported from the #168 word that build.py made
(WORK/cache/pr168.pkl); frames clipped to the future cap and closed forward; donor/recipient pairs = the frozen
data file of #193.
ADAPTED from research/source-assisted-v4/source_aligned_local_v4.py of the outside repository
(github.com/CrocSwap/integer-mult-bounds, Apache-2.0) at commit 187e101.  It was read as text and never run.
op_tags(), local_values() and the body of build() follow op_tags(), local_values() and main() of that script
statement by statement, except its search for pairs and the records it writes; here the pairs are read from a file.
That script came with pull requests #191 and #193 (ikeboy / Avi Eisenberg) and says of itself: "Derived from
research/source-assisted/decision/source_aligned_local.py (icekylinx, GPT-6 Astra assistance, Apache-2.0)."
Changed here: one statement to a line; the graph, the closure and the gauges come from graphgen.py and
compile.py of this directory; the search for pairs and the records written at the end are left out;
load_pairs() is new.  NOTICE, section 6, lists all files of this kind.

usage: python3 -I -B aligned.py pr193      -> WORK/cache/pr193.aligned.pkl
(with the environment variables REBUILD_PAIRS and REBUILD_NAME: that pairs file under that name, see rebuild.py)
"""
import sys, os, json, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f2 import rref, perp, inside
import graphgen, modules, compile as comp
import rbpaths

PR168 = rbpaths.inputs() + "pr168/"
WORK = rbpaths.work()
PAIRS = {"pr193": rbpaths.inputs() + "pr193/research/source-assisted-v4/data/physical-pairs.json"}
if os.environ.get("REBUILD_PAIRS"):       # rebuild.py pairs=<file> name=<name>: another file of hand-over pairs
    PAIRS = {os.environ["REBUILD_NAME"]: os.environ["REBUILD_PAIRS"]}
RB = comp.ROOTBIT


def local_values(g, cut, v):
    rows = [((i, 1),) for i in range(v)]
    for x in range(v, cut):
        aa, bb = g["args"][x]
        d, sign = dict(rows[aa]), g["signs"][x]
        for j, c in rows[bb]:
            z = d.get(j, 0) + sign * c
            if z:
                d[j] = z
            else:
                d.pop(j, None)
        rows.append(tuple(sorted(d.items())))
    return rows


def op_tags(wit):
    args, roots, order = wit["args"], wit["roots"], wit["order"]
    uses = [[] for _ in args]
    for x in order:
        if args[x]:
            for j, y in enumerate(args[x]):
                uses[y].append(2 * x + j)
    for j, r in enumerate(roots):
        uses[r["node"]].append(RB | j)
    incoming = set(wit["arcs"].values())
    tags = []
    for x in order:
        if args[x]:
            tags.append(("sum", x))
        free = [u for u in uses[x] if u not in incoming]
        assert free
        tags.extend(("copy", x, u) for u in free[1:])
    return tags


def load_pairs(path, R, nops, gauged):
    """the hand-over pairs [[donor role, recipient role, operation index or null], ...] of a JSON file, checked for
    form and range against the aligned word: roles below R, operation indices below nops, every recipient a
    gauged role, no role twice.  Whether the hand-overs are LEGAL is not decided here: the frame scan of phys.py
    and gx.check1 at the end of the rebuild decide that."""
    def stop(msg):
        sys.exit("pairs file %s: %s" % (path, msg))
    try:
        data = json.load(open(path))
    except ValueError as e:
        stop("not JSON (%s)" % e)
    pairs = data.get("pairs") if isinstance(data, dict) else None
    if not isinstance(pairs, list):
        stop('no list under the key "pairs"')
    donors, recipients = set(), set()
    for n, q in enumerate(pairs):
        if not (isinstance(q, list) and len(q) == 3 and type(q[0]) is int and type(q[1]) is int
                and (q[2] is None or type(q[2]) is int)):
            stop("entry %d is not [donor, recipient, operation index or null]: %r" % (n, q))
        a, b, t = q
        if not (0 <= a < R and 0 <= b < R):
            stop("entry %d: a role outside 0 .. %d: %r" % (n, R - 1, q))
        if t is not None and not 0 <= t < nops:
            stop("entry %d: an operation index outside 0 .. %d: %r" % (n, nops - 1, q))
        if b not in gauged:
            stop("entry %d: the recipient %d is not a gauged role of the aligned word" % (n, b))
        if a in donors or b in recipients:
            stop("entry %d: its donor or its recipient occurs in an earlier entry: %r" % (n, q))
        donors.add(a)
        recipients.add(b)
    if donors & recipients:
        stop("role %d is a donor in one entry and a recipient in another" % min(donors & recipients))
    return pairs


def build(name):
    t0 = time.time()
    old = pickle.load(open(WORK + "cache/pr168.pkl", "rb"))
    oldg, oldwit, oldword = old["g"], old["wit"], old["word"]
    src = PR168 + "references/paired-cube/sources/"
    local = json.load(open(src + "local_L1.json"))
    local["G"] = {k: "s" for k in local["G"]}
    local["A"] = {k: "fd" for k in local["A"]}
    G = graphgen.Gen(11)
    G.finish(modules.load_triple(src + "tmod_TE_TD_TB3_1_1_4_full_6.0617964e-4.json", 11),
             modules.load_zero_based(src + "pmod_J0_full_6.0666810e-4.json", 45, modules.pair_want(10)),
             modules.load_zero_based(src + "qmod_climb3u_best.json", 9, [511 ^ (1 << i) for i in range(9)]),
             local=local)
    G.fuse("00111100")
    g = G.graph()
    g["kindof"], g["counts"] = G.kindof, G.counts
    h, v, inputs = g["h"], g["v"], g["inputs"]
    oldcut, newcut = v + oldg["counts"]["local"], v + g["counts"]["local"]
    delta = newcut - oldcut
    oldvals, newvals = local_values(oldg, oldcut, v), local_values(g, newcut, v)
    newindex = {row: i + 1 for i, row in enumerate(newvals)}
    node_map = {i + 1: newindex[row] for i, row in enumerate(oldvals) if row in newindex}
    node_map.update({x + 1: x + delta + 1 for x in range(oldcut, len(oldg["args"]))})

    def use_map(u):
        return u if u >> 31 else 2 * node_map[u // 2] + (u & 1)
    newargs, newroots = comp.shift(g)
    arcs = []
    for x, u in json.load(open(PR168 + "references/paired-cube/selected-module/matching-arcs.json")):
        if x <= oldcut or x not in node_map or (not u >> 31 and u // 2 not in node_map):
            continue
        xx, uu = node_map[x], use_map(u)
        val = newroots[uu & (RB - 1)]["node"] if uu >> 31 else newargs[uu // 2][uu & 1]
        if val in newargs[xx]:
            arcs.append([xx, uu])
    prof, wit = comp.closure(g, arcs)
    word = comp.select(g, prof, wit)
    assert all(ca in (-1, 1) and cb in (-1, 1) for ca, cb in word["coef"])
    print("[%s] aligned word: local additions %d (old %d), mapped arcs %d, R=%d, gauges %s  (%.1fs)"
          % (name, g["counts"]["local"], oldg["counts"]["local"], len(arcs), prof["R"], dict(word["gauges"]),
             time.time() - t0), flush=True)
    oldtags, newtags = op_tags(oldwit), op_tags(wit)
    assert len(oldtags) == len(oldword["ops"]) and len(newtags) == len(word["ops"])
    oldframes = [perp(oldwit["ann"][x], h) for _, _, x in oldword["ops"]]
    for i, F in json.load(open(PR168 + "references/paired-cube/physical/frames.json"))["frames"]:
        oldframes[i] = tuple(F)
    desired = {}
    for tag, F in zip(oldtags, oldframes):
        if tag[1] not in node_map:
            continue
        if tag[0] == "copy" and not tag[2] >> 31 and tag[2] // 2 not in node_map:
            continue
        key = ("sum", node_map[tag[1]]) if tag[0] == "sum" else ("copy", node_map[tag[1]], use_map(tag[2]))
        desired[key] = F
    spans = wit["spans"]
    root_frames = {}
    for r, slot in zip(newroots, word["rootroles"]):
        root_frames[slot] = spans[r["node"]] if r["kind"] == "center" else \
            perp(rref(inputs[t] for t in r["targets"]), h)
    R, ops, pset = prof["R"], word["ops"], set(word["phase1"])
    previous, parents = [-1] * R, []
    for i, (aa, bb, x) in enumerate(ops):
        parents.append((previous[aa], previous[bb]))
        previous[aa] = previous[bb] = i
        if x <= newcut and len(spans[x]) <= 2:
            pset.add(i)
    todo = list(pset)
    while todo:
        i = todo.pop()
        for j in parents[i]:
            if j >= 0 and j not in pset:
                pset.add(j)
                todo.append(j)
    chron = sorted(pset) + [i for i in range(len(ops)) if i not in pset]
    cap = [perp(root_frames[s], h) if s in root_frames else () for s in range(R)]
    maxima = [None] * len(ops)
    for i in reversed(chron):
        aa, bb, x = ops[i]
        A = rref(cap[aa] + cap[bb])
        maxima[i] = perp(A, h)
        cap[aa] = cap[bb] = A
        assert inside(spans[x], maxima[i])
    selected = [z for z in word["selected"] if z["rank"] == h - 4]
    gauge = {z["role"]: perp(tuple(z["annihilator"]), h) for z in selected}
    current = [() for _ in range(R)]
    for x, slot in word["sources"].items():
        current[slot] = (inputs[int(x) - 1],)
    for slot, U in gauge.items():
        current[slot] = U
    frames, mapped, minimal = [None] * len(ops), 0, 0
    for i in chron:
        aa, bb, x = ops[i]
        if x <= newcut:
            want = spans[x]
        else:
            want = desired.get(newtags[i], spans[x])
            mapped += newtags[i] in desired
            want = perp(rref(perp(want, h) + perp(maxima[i], h)), h)
        U = rref(spans[x] + want + current[aa] + current[bb])
        assert inside(U, maxima[i]), ("forward frame outside future cap", i)
        frames[i] = current[aa] = current[bb] = U
        if x <= newcut and len(spans[x]) == 2:
            assert U == spans[x]
            minimal += 1
    for slot, U in root_frames.items():
        assert inside(current[slot], U)
    pairs = load_pairs(PAIRS[name], R, len(ops), set(gauge))
    kept = {b for _, b, _ in pairs}
    word["selected"] = [z for z in selected if z["role"] in kept]
    word["phase1"] = sorted(pset)
    moved = [[i, list(U)] for i, U in enumerate(frames) if U != perp(wit["ann"][ops[i][2]], h)]
    print("    phase one %d ops; mapped query ops %d; minimal local pair ops %d; moved frames %d; pairs %d; gauges "
          "kept %d of %d  (%.1fs)" % (len(pset), mapped, minimal, len(moved), len(pairs), len(word["selected"]),
                                      len(selected), time.time() - t0), flush=True)
    pickle.dump(dict(g=g, prof=prof, wit=wit, word=word, frames=frames, moved=moved, pairs=pairs,
                     root_frames=root_frames, newcut=newcut),
                open(WORK + "cache/%s.aligned.pkl" % name, "wb"))


if __name__ == "__main__":
    for nm in sys.argv[1:]:
        build(nm)
