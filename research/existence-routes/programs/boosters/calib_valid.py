# e3-boost calib_valid.py: every new certificate (B, C_k, distance sets, projections) must pass relations that are
# known to be realisable on registers of m <= 8 rows: plain [k] <= k[1]; Kronecker [a+b] <= [a] + [b];
# undo [k] <= [k'] + [k' - k] (H on k' rows, then H again on the k' - k extra rows); A5 [1] <= 2[m'] (m' <= m).
import certlib as C
rels = []
for m in range(2, 9):
    for k in range(2, m + 1): rels.append((f"plain [{k}]<={k}[1]", m, C.mk({k: 1}, {1: k})))
    for a in range(1, m):
        for b in range(a, m - a + 1): rels.append((f"Kron [{a+b}]<=[{a}]+[{b}]", m, C.mk({a + b: 1}, {a: 1, b: 1} if a != b else {a: 2})))
    for k in range(1, m):
        for kk in range(k + 1, m + 1): rels.append((f"undo [{k}]<=[{kk}]+[{kk-k}]", m, C.mk({k: 1}, {kk: 1, kk - k: 1} if kk - k != kk else {kk: 2})))
    for mm in range(2, m + 1): rels.append((f"A5 [1]<=2[{mm}]", m, C.mk({1: 1}, {mm: 2})))
bad = 0; tests = 0
for name, m, v in rels:
    res = C.certify(v, m); tests += len(res)
    fails = [r for r in res if r[3] < 0]
    if fails: bad += 1; print("FAIL", name, "register", m, [(n, b, s) for n, b, s, _ in fails][:4])
print(f"{len(rels)} realisable relations on registers 2..8, {tests} certificate tests, relations failing some test: {bad}")
