"""prune.py <hierarchy.json> <drop.json>: compile the hierarchy and drop the target stops that do not pay."""
import sys, json
sys.dont_write_bytecode = True
from ilib import *
h = load(sys.argv[1])
c, st, drop, w = prune_low(h)
f, q = figure(c)
print("compiler + pruning: R %d low %d figure %d" % (st["R"], st["lowreads"], f))
json.dump([list(x) for x in sorted(drop)], open(sys.argv[2], "w"))
