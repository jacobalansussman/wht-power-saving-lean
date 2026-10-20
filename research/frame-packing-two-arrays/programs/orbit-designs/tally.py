#!/usr/bin/env python3
"""tally.py (pk2-orbits): one markdown table over all orb.py outputs given, grouped by target (mode, block) and run
family (orbit kinds, char), with summ.py's classification of each design's final verdicts."""
import sys, re, ast, collections, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summ import classify

TARGET = {("two", 4): "2 in 7, block 4", ("two", 5): "2 in 9, block 5", ("two2", 5): "2 in 8, block 5",
          ("three", 4): "3 in 10, block 4", ("two", 3): "control: 2 in 5, block 3", ("two2", 4): "control: 2 in 6, block 4"}
COLS = ["FOUND", "one output exact", "dead for every scalar choice (relaxations)",
        "fewer than W states reach level m (dead for every scalar)", "exact no at drawn scalars (full solve)",
        "exact no at drawn scalars (layered screen)", "exact no at drawn scalars (psi at those scalars)", "cap (undecided)"]
groups = collections.OrderedDict()
for fn in sys.argv[1:]:
    f = open(fn); hdr = f.readline()
    m = re.match(r"orb.py mode=(\S+) m=(\d+) .*? kinds=(\[.*\])", hdr)
    if not m: continue
    mode, blk, kinds = m.group(1), int(m.group(2)), ast.literal_eval(m.group(3))
    fam = "all kinds" if len(kinds) == 5 else ",".join(kinds)
    if "char" in fn: fam += " (char scalars)"
    wave = os.path.basename(fn).split("_")[1] if os.path.basename(fn).startswith("out_") else "?"
    wave = wave if wave in ("B", "C") else "A"
    key = (TARGET.get((mode, blk), f"{mode} {blk}"), fam)
    g = groups.setdefault(key, collections.Counter())
    for line in f:
        mm = re.match(r"trial (\d+): a=(.*?) dup=(\d+) orbits=(\[.*?\]) hist=(\[.*\])$", line.strip())
        if not mm: continue
        orbits = ast.literal_eval(mm.group(4)); h = ast.literal_eval(mm.group(5))[-1]
        g["designs"] += 1; g["orbits"] += len(orbits); g[classify(h)] += 1; g["wave " + wave] += 1
short = ["found", "one exact", "dead, every scalar", "< W top states", "no (full solve)", "no (layered)",
         "no (psi, drawn scalars)", "cap"]
print("| target | orbit kinds | designs (waves A/B/C) | orbits | " + " | ".join(short) + " |")
print("|---|---|---|---|" + "---|" * len(short))
tot = collections.Counter()
for (tgt, fam), g in sorted(groups.items(), key=lambda kv: (kv[0][0].startswith("control"), kv[0])):
    waves = "/".join(str(g.get("wave " + w, 0)) for w in "ABC")
    print(f"| {tgt} | {fam} | {g['designs']} ({waves}) | {g['orbits']} | " + " | ".join(str(g.get(c, 0)) for c in COLS) + " |")
    if not tgt.startswith("control") and fam != "cycle": tot.update(g)
print("\n(the pure 'cycle' family places no orbit at all: before a first landing each input has one junk term, so no"
      " cycle exists; it is left out of the totals)")
print(f"\ntargets total: {tot['designs']} designs, {tot['orbits']} orbit injections; " +
      ", ".join(f"{s}: {tot.get(c, 0)}" for s, c in zip(short, COLS)))
