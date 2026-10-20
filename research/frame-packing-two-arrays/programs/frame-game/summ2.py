"""summ2.py: per file, designs; dead for every scalar choice (linearised psi refuses some output); exact 'no' for the
drawn scalars (some output inconsistent or unproduced); undecided (every output yes or cap, and linearised psi passes
all outputs); one output exact; unit candidates."""
import sys, ast, collections
tot = collections.Counter()
for f in sys.argv[1:]:
    txt = open(f).read().splitlines()
    c = collections.Counter()
    for l in txt:
        l = l.strip()
        if not (l.startswith("(") and l[-1].isdigit()): continue
        try:
            k, v = l.rsplit(" ", 1); k = ast.literal_eval(k); v = int(v)
        except Exception: continue
        st = list(k[-1])
        pl = [x for x in st if str(x).startswith("psilin")]
        st = [x for x in st if not str(x).startswith("psilin")]
        c["designs"] += v
        dead_all = bool(pl) and "0" in pl[0].split(":")[1]
        dead_drawn = any(("linear" in str(s)) or ("target" in str(s)) for s in st)
        if dead_all: c["dead for all scalars"] += v
        if dead_drawn: c["exact no (drawn scalars)"] += v
        if not dead_all and not dead_drawn: c["undecided (cap) "] += v
        if any(s == "yes" for s in st): c["one output exact"] += v
        if all(s == "yes" for s in st): c["BOTH"] += v
    c["candidates"] = sum(1 for l in txt if "CANDIDATE" in l)
    head = next((l for l in txt if l.startswith("design")), f)
    print(head); print("   ", dict(c))
    tot.update(c)
print("TOTAL", dict(tot))
