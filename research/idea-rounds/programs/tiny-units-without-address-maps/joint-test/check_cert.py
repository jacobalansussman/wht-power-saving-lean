#!/usr/bin/python3
"""joint/proof/check_cert.py -- checks static_cert_W3_k8.json from the label formulas alone (no flow, no lib.py)."""
import json, itertools
def dot(a, b): return bin(a & b).count("1") & 1
labels = [(x, z) for x in range(8) for z in range(8)]
blk = lambda x, z: (x, dot(x, z)) if x else (0, 0)             # block of the label (x|z)
dirl = lambda w, x, z: ((x ^ w) if dot(w, z) == 0 else x, z)   # dir(w) on labels
kerl = lambda x, z: (x ^ 7 ^ z, z)                             # the kernel on labels
seen = set()
import os
CERT = os.environ.get("CERT_JSON", "static_cert_W3_k8.json")
if not os.path.exists(CERT):   # default: the copy in this repository's data folder
    CERT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "data",
                        "tiny-units-without-address-maps", "joint-test", "static_cert_W3_k8.json")
for ms, B, lhs, rhs in json.load(open(CERT)):
    B = set(map(tuple, B))
    need = 3 * sum(1 for (x, z) in labels if blk(x, z) not in B and blk(*kerl(x, z)) in B)
    got = sum(sum(1 for (x, z) in labels if blk(x, z) not in B and blk(*dirl(w, x, z)) in B) for w in ms)
    assert got == lhs and need == rhs and got < need and len(ms) == 8
    seen.add(tuple(ms))
assert seen == set(itertools.combinations_with_replacement(range(1, 8), 8))
print("static_cert_W3_k8.json: all %d multisets of 8 directions have a set B with (labels the moves can carry in) < 3 x (labels the kernel carries in): CHECKED" % len(seen))
