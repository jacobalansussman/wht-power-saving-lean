"""split_tight.py: split a survivor dump into classes with and without a tight form (|D_Q| = r_Q >= 1)."""
import sys
from xbase import Forms, System
from tile_alt import target_matrices
m = int(sys.argv[1]); path = sys.argv[2]
mats, nonalt = target_matrices(m); fo = Forms(m)
a = open(path + ".notight", "w"); b = open(path + ".tight", "w"); na = nb = 0
for L in open(path):
    if "|" not in L: continue
    d, t = L.strip().split("|"); counts = {}
    for z in d.split(","): counts[int(z)] = counts.get(int(z), 0) + 1
    sy = System(fo, [mats[int(x)] for x in t.split(",")], counts)
    if len(sy.tight()): b.write(L); nb += 1
    else: a.write(L); na += 1
print(f"{path}: {na} classes with no tight form, {nb} with one")
