"""an-synth calc/twosided.py -- what a second unpaid coordinate (two-sided slack) would cost inside our U8 local word.
Our own code, stdlib, reads our own lowering (k-lower/out/U8X-parsed.json.gz, gcert/1, format as in tools/gx/gxcore.py).
Model (argued in NOTE section 4): every target y_S also starts on a line v_S, with v-lines a mirrored copy of the port
lines (v_S . v_T = u_S . u_T over F_2, in coordinates disjoint from the u-window), and every source x_T stays inside
v_T^perp. A register's frame only grows, so once a slot has met target y_S its frame contains v_S; a later gate with
source x_T then needs v_S . v_T = 0. Each violation forces the slot to drop v-span (a descent plus the re-climb).
This counts, per slot, the conflicting v's (y met before an x it is not orthogonal to) and the F_2 rank of their span:
a first-order estimate of the extra paid rank per call. It also counts the scatter: each retained total is read by the
targets of its scatter column; under two-sided labels the total's frame must contain their v's afterwards.
Run: python3 -I calc/twosided.py"""
import gzip, json, os, sys
P = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), "../../data/inputs/U8X-parsed.json.gz"))
c = json.load(gzip.open(P, "rt"))
v, R, ports = c["v"], c["R"], c["ports"]
def cls(r): return "x" if r < v else ("y" if r < 2 * v else "s")
def regs(g): return [g[2]] + [t for t, _, _ in g[3]] + (list(g[4]) if g[0] == "out" else [])
def par(x): return bin(x).count("1") & 1
def rank_f2(vs):
    basis = {}
    for r in vs:
        while r:
            h = r.bit_length() - 1
            if h in basis: r ^= basis[h]
            else: basis[h] = r; break
    return len(basis)
gates = [("A", g) for g in c["A"]] + [("B", g) for g in c["B"]]
# scatter: after A, retained totals (slots) are read into targets per table row t -> list of (k, num, den)
ret_slot = {k: s for k, s, f in c["ret"]}
scat_rows = c["scat"]["table"]
seq = {}  # slot -> list of events ('x', port) / ('y', port)
def push(r, ev): seq.setdefault(r, []).append(ev)
for ph, g in gates:
    rs = regs(g)
    xs = [r for r in rs if cls(r) == "x"]
    ys = [r for r in rs if cls(r) == "y"]
    for r in rs:
        if cls(r) != "s": continue
        for x in xs: push(r, ("x", x))
        for y in ys: push(r, ("y", y - v))
    if ph == "A" and g is c["A"][-1]:
        pass
# insert the scatter between A and B: rebuild with an explicit marker
seq = {}
for g in c["A"]:
    rs = regs(g)
    for r in rs:
        if cls(r) == "s":
            for x in rs:
                if cls(x) == "x": push(r, ("x", x))
            for y in rs:
                if cls(y) == "y": push(r, ("y", y - v))
for t, row in enumerate(scat_rows):
    for k, a, b in row:
        push(ret_slot[k], ("y", t))
for g in c["B"]:
    rs = regs(g)
    for r in rs:
        if cls(r) == "s":
            for x in rs:
                if cls(x) == "x": push(r, ("x", x))
            for y in rs:
                if cls(y) == "y": push(r, ("y", y - v))
tot_slots = len(seq)
conf_slots, extra_rank, total_conf_pairs = 0, 0, 0
totals_extra = 0
for r, ev in seq.items():
    seen_y = set(); bad_v = set()
    for kind, p in ev:
        if kind == "y":
            seen_y.add(p)
        else:
            for q in seen_y:
                if par(ports[q] & ports[p]):   # v_q . v_p = u_q . u_p = 1: conflict
                    bad_v.add(q); total_conf_pairs += 1
    if bad_v:
        conf_slots += 1
        rk = rank_f2([ports[q] for q in bad_v])
        extra_rank += 2 * rk
        if r in ret_slot.values(): totals_extra += 2 * rk
# direct x -> y gates (x and y at one frame) need u_x . u_y = 0 (already) and v_x . v_y = 0 (same under mirroring)
direct_bad = 0
for ph, g in gates:
    rs = regs(g)
    xs = [r for r in rs if cls(r) == "x"]; ys = [r - v for r in rs if cls(r) == "y"]
    direct_bad += sum(1 for x in xs for y in ys if par(ports[x] & ports[y]))
out = dict(slots_with_events=tot_slots, slots_needing_a_v_drop=conf_slots, conflict_events=total_conf_pairs,
           extra_paid_rank_per_call_estimate=extra_rank, of_which_retained_totals=totals_extra,
           direct_xy_gates_violating=direct_bad)
print(json.dumps(out, indent=1))
w, cst = 99.5556, 64.0
D2 = 2 * w - cst - extra_rank * (w / 224 * 2) / 2  # rescale: extra rank per call is per local call; D is per call too
print("per call (template average): one-sided D = %.2f ; two-sided D ~ 2w - cst - extra = %.2f - %.0f - %d = %.2f" %
      (w - cst, 2 * w, cst, extra_rank, 2 * w - cst - extra_rank))
json.dump(out, open(os.path.join(os.path.dirname(__file__), "twosided.json"), "w"), indent=1)
