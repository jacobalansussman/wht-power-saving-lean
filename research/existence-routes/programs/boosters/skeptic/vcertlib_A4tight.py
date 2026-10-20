# v-e3-boost: which of e3-boost's own certlib tests kill the block-4 A-tight relation [1] + 5[3] <= 4[2] + 2[4]?
import sys
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
import certlib as C
v = {1: -1, 2: 4, 3: -5, 4: 2}
ver, bad, res = C.verdict(v, 4)
kinds = {}
for n, b, s, x in bad: kinds.setdefault((n.rstrip("0123456789S") if n.startswith("Dist") else n, s), 0); kinds[(n.rstrip("0123456789S") if n.startswith("Dist") else n, s)] += 1
print(ver, len(res), "tests;", len(bad), "failing; by family/status:", kinds)
print("e.g.", [(n, b, s, float(x)) for n, b, s, x in bad if n.startswith("Dist4")][:6], [(n, s, float(x)) for n, b, s, x in bad if n.startswith("hook")])
