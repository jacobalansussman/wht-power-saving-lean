#!/usr/bin/env python3
"""growth_summ.py: per repair rule, mean orphan count and mean number of terms after r repairs, and how many
designs had hit the closure cap ('c') by repair r (growth.py outputs)."""
import sys, re, collections
for fn in sys.argv[1:]:
    data = {"single": [], "orbit": []}
    for line in open(fn):
        m = re.match(r"trial \d+ a=\S+ (single|orbit)\s*: (.*?)(  orbits=.*)?$", line.rstrip())
        if not m: continue
        seq = [tuple(x.split("/")) for x in m.group(2).split()]
        data[m.group(1)].append(seq)
    print(f"== {fn}")
    for rule, seqs in data.items():
        out = []
        for r in range(7):
            alive = [s for s in seqs if len(s) > r]
            if not alive: break
            orph = [int(s[r][1]) for s in alive if s[r][1] != "-"]
            terms = [int(s[r][2]) for s in alive]
            capd = sum(1 for s in seqs if any("c" in x[0] for x in s[:r + 1]))
            out.append(f"r={r}: n={len(alive)} orphans {sum(orph)/len(orph) if orph else float('nan'):.1f} "
                       f"terms {sum(terms)/len(terms):.0f} capped-by-now {capd}/{len(seqs)}")
        print(f"   {rule}:"); [print("      " + o) for o in out]
