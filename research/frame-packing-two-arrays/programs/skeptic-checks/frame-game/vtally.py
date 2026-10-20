"""vtally.py (v-fp-build): my own recount of fp-build's production tallies from the raw Colab output lines.
A line is '<key tuple> <count>' whose last tuple element is the statuses (one per output, then 'psilin:bits')."""
import re, sys, ast, collections
import os
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "data", "frame-game") + os.sep
files = [("2 in 7, block 4", "design-runs/out_design_4_1500_12_61_psilin.txt"),
         ("2 in 7 dup, block 4", "design-runs/out_design_4_1500_12_62_dup_psilin.txt"),
         ("2 in 9, block 5", "design-runs/out_design_5_600_10_71_psilin.txt"),
         ("2 in 8, block 5", "design-runs/out_design_5_600_10_72_t2_psilin.txt"),
         ("2 in 11, block 6", "design-runs/out_design_6_300_8_81_psilin.txt"),
         ("2 in 10, block 6", "design-runs/out_design_6_300_8_82_t2_psilin.txt"),
         ("3 in 10, block 4", "design-runs/out_designw_3in10_b4_s91.txt"),
         ("ctl 3 in 7, block 3", "design-runs/out_designw_3in7_b3_s92.txt"),
         ("chk 2 in 3, block 2", "out_design_m2_s3.txt")]
pat = re.compile(r"^\s*(\(.*\))\s+(\d+)\s*$")
T = collections.Counter()
print(f"{'target':22s} {'designs':>7s} {'both':>5s} {'deadAll':>7s} {'noDrawn':>7s} {'open':>5s} {'oneYes':>6s} {'cand':>4s}")
for name, f in files:
    c = collections.Counter(); cand = 0
    for line in open(base + f):
        if "CANDIDATE" in line: cand += 1
        mt = pat.match(line)
        if not mt: continue
        key = ast.literal_eval(mt.group(1)); v = int(mt.group(2))
        st = key[-1]
        outs = [x for x in st if not x.startswith("psilin")]
        bits = [x for x in st if x.startswith("psilin")][0].split(":")[1]
        c["designs"] += v
        dead_all = "0" in bits
        drawn_no = any(x.startswith("no") for x in outs)
        if all(x == "yes" for x in outs): c["both"] += v
        if dead_all: c["deadAll"] += v
        if drawn_no: c["noDrawn"] += v
        if not dead_all and not drawn_no: c["open"] += v
        if any(x == "yes" for x in outs): c["oneYes"] += v
    print(f"{name:22s} {c['designs']:7d} {c['both']:5d} {c['deadAll']:7d} {c['noDrawn']:7d} {c['open']:5d} {c['oneYes']:6d} {cand:4d}")
    if not name.startswith(("ctl", "chk")): T.update(c)
print("targets total:", dict(T), f" deadAll share {T['deadAll']/T['designs']:.3f}, open share {T['open']/T['designs']:.3f}")
