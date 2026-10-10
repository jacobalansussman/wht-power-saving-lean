#!/usr/bin/python3 -I
"""phys.py: step 4 of the rebuild.  Applies the layer (moved frames, reuse pairs) to the word that rewrite.py
made, builds the literal stage word, scans the frames of every register, counts the blocks and the unit moves,
and classifies every block.

usage: python3 -I -B phys.py pr193xword        (after rewrite.py)
writes WORK/cache/<name>.stage.pkl (the layer, the block counts and N; emit.py reads it) and
WORK/logs/phys-<name>.json
"""
import sys, os, json, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter
from f2 import kind, step_kind
import stage
import rbpaths

WORK = rbpaths.work()


def layer_of(name):
    """the layer that rewrite.py wrote for this word"""
    return pickle.load(open(WORK + "cache/%s.layer.pkl" % name, "rb"))["layer"]


def summary(g, prof, st, H, blocks):
    h, v = g["h"], g["v"]
    Rp = len(st["live"])
    tails = Counter(len(st["gauge"][s]) for s in st["live"] if s in st["gauge"])
    C = Counter()
    for c in "xysc":
        for r, n in H[c].items():
            C[r] += 3 * n
    C.update({3 * d: n for d, n in tails.items() if d})
    C[2] += 2 * v
    m, W = 3 * h, 2 * v + Rp
    mass = sum(r * n for r, n in C.items())
    N = sum(len(B) - len(A) for _, A, B, _ in blocks)
    return dict(h=h, v=v, R_word=prof["R"], pairs=len(st["pairs"]), sinks=len(st["sink"]), physical_R=Rp,
                W_per_vertex=W, m=m, rank_per_vertex=mass, deficit_per_vertex=W * m - mass,
                deficit_formula=2 * v - 3 * prof["loss"], loss=prof["loss"], N_moves=N, events=len(st["events"]),
                moved_frames=st["moved"], gauge_tails={str(k): n for k, n in tails.items()},
                H={c: {str(r): n for r, n in sorted(H[c].items())} for c in "xysc"},
                child_histogram={str(r): n for r, n in sorted(C.items())})


def classify(g, st, blocks):
    """every block: class of register x / y / s (slot) / c (star copy), what closes it, and whether the step is
    legal in the label calculus of this repository"""
    h, v, ev = g["h"], g["v"], st["events"]
    kindof = g["kindof"]
    memo = {}
    tab = Counter()
    mass = Counter()
    fk = {}
    for r, A, B, n in blocks:
        cls = "c" if r < 0 else "x" if r < v else "y" if r < 2 * v else "s"
        if n >= len(ev):
            top = "final"
        else:
            e = ev[n]
            top = e["k"]
            if top == "op":
                top = "op:" + kindof.get(e["node"] - 1, "port")
            elif top in ("undo", "uninj"):
                top = "final"
        if (A, B) not in memo:
            memo[A, B] = step_kind(A, B, h)
        for F in (A, B):
            if F not in fk:
                fk[F] = kind(F, h)
        tab[cls, top, memo[A, B]] += 1
        mass[cls, memo[A, B]] += len(B) - len(A)
    frames = Counter(fk.values())
    return tab, mass, frames


def main(name):
    t0 = time.time()
    d = pickle.load(open(WORK + "cache/%s.pkl" % name, "rb"))
    g, prof, wit, word = d["g"], d["prof"], d["wit"], d["word"]
    layer = layer_of(name)
    out = {}
    for tag, us in (("no sinks", False), ("with sinks", True)):
        if us and not layer.get("sinks"):
            continue
        st = stage.build(g, wit, word, layer, use_sinks=us)
        blocks, H = stage.framescan(st)
        S = summary(g, prof, st, H, blocks)
        out[tag] = S
        print("[%s, %s] R_word=%d pairs=%d sinks=%d physical R=%d  W=%d  rank=%d  deficit=%d (2v-3loss=%d)  N=%d  "
              "events=%d moved frames=%d blocks=%d  (%.1fs)"
              % (name, tag, S["R_word"], S["pairs"], S["sinks"], S["physical_R"], S["W_per_vertex"],
                 S["rank_per_vertex"], S["deficit_per_vertex"], S["deficit_formula"], S["N_moves"], S["events"],
                 S["moved_frames"], len(blocks), time.time() - t0), flush=True)
    last = "with sinks" if "with sinks" in out else "no sinks"
    st = stage.build(g, wit, word, layer, use_sinks=(last == "with sinks"))
    blocks, H = stage.framescan(st)
    tab, mass, frames = classify(g, st, blocks)
    print("    distinct frames at block ends: %s" % dict(frames))
    tot = Counter()
    for (cls, top, lg), n in tab.items():
        tot[cls, lg] += n
    print("    blocks by class and legality: %s" % {"%s/%s" % k: n for k, n in sorted(tot.items())})
    print("    moves by class and legality:  %s" % {"%s/%s" % k: n for k, n in sorted(mass.items())})
    for (cls, top, lg), n in sorted(tab.items()):
        print("      %s %-12s %-5s %7d" % (cls, top, lg, n))
    out["classification"] = {"%s|%s|%s" % k: n for k, n in sorted(tab.items())}
    out["moves_by_class_legality"] = {"%s|%s" % k: n for k, n in sorted(mass.items())}
    out["frames_at_block_ends"] = dict(frames)
    json.dump(out, open(WORK + "logs/phys-%s.json" % name, "w"), indent=1, sort_keys=True)
    pickle.dump(dict(layer=layer, summary=out), open(WORK + "cache/%s.stage.pkl" % name, "wb"))
    print("    done (%.1fs)" % (time.time() - t0), flush=True)


def stop(name, e):
    """a check of the stage word (an assert of stage.py) failed: one plain line and exit code 1, in place of a
    Python traceback.  With rebuild.py pairs=<file> this is where a pairs file of the right form ends when its
    hand-overs are not legal."""
    import traceback
    at = traceback.extract_tb(e.__traceback__)[-1]
    why = e.args[0] if e.args else at.line
    if isinstance(why, tuple) and len(why) == 5 and why[0] == "frame step not nested":
        why = ("frame step not nested: an event of kind %r needs register %d at a frame of dimension %d, and the "
               "register stands at a frame of dimension %d that is not inside it" % (why[1], why[2], why[4], why[3]))
    msg = ("phys.py %s: a check of the stage word failed (%s, line %d): %s.  The rebuild stops here; no "
           "certificate is written." % (name, os.path.basename(at.filename), at.lineno, why))
    if os.environ.get("REBUILD_PAIRS"):
        msg += ("  The pairs file %s has the right form, but its hand-overs are not legal for this word."
                % os.environ["REBUILD_PAIRS"])
    sys.exit(msg)


if __name__ == "__main__":
    for nm in sys.argv[1:]:
        try:
            main(nm)
        except AssertionError as e:
            stop(nm, e)
