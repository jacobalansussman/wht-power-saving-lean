# v-e3-boost vvalid.py: my certificate families against relations known to be realisable, on registers M >= m (b = 1).
import vsurv as S
rels = {}
for k in range(2, 8): rels[f"plain[{k}]<={k}[1]"] = {k: -1, 1: k}
for m in range(2, 8): rels[f"A5 [1]<=2[{m}]"] = {1: -1, m: 2}
for a, b in [(1, 1), (1, 2), (2, 2), (2, 3), (3, 4)]:
    v = {a + b: -1}; v[a] = v.get(a, 0) + 1; v[b] = v.get(b, 0) + 1; rels[f"Kron[{a+b}]<=[{a}]+[{b}]"] = v
for k in range(1, 6):
    for kp in range(k + 1, 8):
        v = {k: -1}; v[kp] = v.get(kp, 0) + 1; v[kp - k] = v.get(kp - k, 0) + 1
        rels[f"undo[{k}]<=[{kp}]+[{kp-k}]"] = {d: c for d, c in v.items() if c}
fails = {}
for name, v in rels.items():
    m = max(v)
    for M in range(m, min(m + 4, 9) + 1):
        res = S.tests(v, M, Mextra=0)
        for n, x in res:
            if x < 0 and "x2" not in n: fails.setdefault(n.split("S")[0].split("T")[0] if n.startswith(("Dist", "Lines")) else n, []).append((name, n, int(x)))
print(len(rels), "relations; certificate rows failing on a realisable relation:")
for k, l in fails.items(): print(" ", k, len(l), l[:3])
