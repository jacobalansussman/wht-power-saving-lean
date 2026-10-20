# lockstep_sim.py (meta/ex-laser, 14 Oct 2026).  Literal simulation of the lockstep schedule of PROOF.md Lemma 2
# (no deadlock): every unit issues its M paid moves in order, and move k may be issued only after every job created
# by its moves 1..k-1 has been completed (i.e. the units that consumed those jobs have finished).  Exact fractions.
# Classes: a class is a set of units at one generation g = a+b (column c = 3^a 2^b) with scale lam; phase k of a class
# issues move k of all its units; the resulting fluid, plus lam*Source/M of frontier fluid, feeds the child class at
# generation g-1 with scale lam/M.  The simulation checks at every class: (i) the fluid it receives covers exactly
# the inputs of its units (no deficit, no surplus except size-1 leaves), (ii) all inputs were issued before the class
# starts (clock), (iii) totals: sources used = Source, leaves = M*u(0,0) (times the normalisation).
import sys
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from fractions import Fraction as Fr
from math import comb
from frontier_flow import flow

def run(L, W3, W2, M):
    U, src, leaves, bad = flow(L, W3, W2, M, exact=True)
    assert not bad
    size = lambda a, b: 3 ** a * 2 ** b
    pt = {}                                   # size -> (a, b)
    for a in range(0, 40):
        for b in range(0, 60):
            if size(a, b) <= L: pt[size(a, b)] = (a, b)
    gen = lambda s: sum(pt[s])
    gmax = max(a + b for (a, b) in U)
    units_at = lambda g: {ab: u for ab, u in U.items() if sum(ab) == g}
    src_at = lambda g: {s: w for s, w in src.items() if gen(s) == g}
    clock = [0]; stats = dict(classes=0, src_used=Fr(0), leaves=Fr(0), max_def=Fr(0))
    issued = {}                                # (class id, k) -> time issued

    def process(g, lam, fluid, cid):
        # fluid: dict size -> amount available at generation g+1 (produced by the parent phase + its sources)
        stats['classes'] += 1
        start = clock[0]; clock[0] += 1
        need = {}
        for (a, b), u in units_at(g).items():
            c = size(a, b)
            need[3 * c] = need.get(3 * c, 0) + W3 * lam * u
            need[2 * c] = need.get(2 * c, 0) + W2 * lam * u
        for s in set(need) | set(fluid):
            d = fluid.get(s, 0) - need.get(s, 0)
            if d != 0:
                raise AssertionError(f"class {cid} gen {g}: size {s} surplus/deficit {d}")
        for k in range(1, M + 1):
            # issue move k of every unit in this class: fluid at generation g
            issued[(cid, k)] = clock[0]; clock[0] += 1
            out = {}
            for (a, b), u in units_at(g).items():
                out[size(a, b)] = out.get(size(a, b), 0) + lam * u
            add = {s: lam * w / M for s, w in src_at(g).items()}
            stats['src_used'] += sum(add.values())
            for s, w in add.items(): out[s] = out.get(s, 0) + w
            if 1 in out:                       # size-1 jobs are leaves (one plain pass each)
                stats['leaves'] += out.pop(1)
            if g >= 1 and out:
                process(g - 1, lam / M, out, cid + (k,))
            elif out:
                raise AssertionError("fluid left below generation 0")
            # all jobs of move k are now complete; the next move may be issued
        return

    # root: units at gmax are fed only by sources at generation gmax+1
    top = {s: w for s, w in src.items() if gen(s) == gmax + 1}
    stats['src_used'] += sum(top.values())
    process(gmax, Fr(1), top, ())
    return stats, sum(src.values()), leaves

if __name__ == "__main__":
    for (W3, W2, M, L) in ((1, 4, 10, 12), (1, 4, 10, 30), (1, 1, 4, 30), (1, 4, 10, 54), (7, 38, 88, 30)):
        st, S, lv = run(L, W3, W2, M)
        print(f"{W3}+{W2} in {M}, L={L}: classes {st['classes']}, every class exactly fed (no surplus/deficit), "
              f"sources used {st['src_used']} = Source total {S}: {st['src_used'] == S}; "
              f"leaves {st['leaves']} = M*u(0,0) {lv}: {st['leaves'] == lv}", flush=True)
