"""v-xc3-s1s2-compute extra: (1) chi S_last(row 1) chi^-1 vs P; (2) E^-1 Dcz E Dcz^-1; (3) Theorem S at lambda = 1
for letters regular at 1 (q = 1): symbols over F2, rank growth of words with s(w) != I; k, l, P, chi must NOT be
regular at 1."""
import sys, random
sys.argv = ["x"]
src = open("vcheck.py").read().split('print("== (A) regularity')[0]
exec(src)
random.seed(7)
CACHE = {}
def letter_M(nm, f):
    key = (nm, f)
    if key not in CACHE:
        base = nm.replace("^-1", "")
        X = LETTERS[base][1](f)
        if nm.endswith("^-1"): X = f2inv(X)
        CACHE[key] = X
    return CACHE[key]
def word_M(w, f):
    M = np.eye(4 * f, dtype=np.uint8)
    for nm in w: M = f2mul(M, letter_M(nm, f))
    return M
def rk(w, f): return f2rank(word_M(w, f) ^ np.eye(4 * f, dtype=np.uint8))
LETTERS["Slast1"] = ("free", lambda f: label_from_phase(full_B(f, lambda f: [(bidx(f - 1, 1), bidx(f - 1, 1))]), f))
for w in (["chi", "Slast1", "chi^-1", "P^-1"], ["chi_row1", "Slast1", "chi_row1^-1", "P^-1"]):
    print("rank(M(w) - I), f = 12, 20, 33:", [rk(w, f) for f in (12, 20, 33)], " ".join(w))
w = ["E^-1", "Dcz", "E", "Dcz^-1"]
f = 20; M = word_M(w, f)
print("E^-1 Dcz E Dcz^-1: is phase type (x-part identity, z += B x):",
      np.array_equal(M[np.ix_([4*j+t for j in range(f) for t in (0,1)], [4*j+t for j in range(f) for t in range(4)])],
                     np.eye(4*f, dtype=np.uint8)[[4*j+t for j in range(f) for t in (0,1)]]))
Bsub = M[np.ix_([4*j+t for j in range(f) for t in (2,3)], [4*j+t for j in range(f) for t in (0,1)])]
nzs = sorted({(a // 2 - b // 2, a % 2, b % 2) for a, b in zip(*np.nonzero(Bsub)) if 3 < a // 2 < f - 4})
print("   z += B x couplings in the bulk (column offset, z-row, x-row):", nzs)
# lambda = 1
one_pow = lambda d: 1
ok = True; SYM1 = {}
for name, (kind, fn) in LETTERS.items():
    if kind == "control": continue
    ns = [deep_toeplitz(fn(f), {0: 1}, f) for f in (23, 31, 40)]
    nsi = [deep_toeplitz(f2inv(fn(f)), {0: 1}, f) for f in (23, 31, 40)]
    reg = all(n is not None for n in ns + nsi) and all(set(ns[0]) == set(x) and all(np.array_equal(ns[0][d], x[d]) for d in ns[0]) for x in ns)
    if reg:
        s = symbol_at(ns[0], {0: 1}, one_pow); si = symbol_at(nsi[0], {0: 1}, one_pow)
        SYM1[name] = s; SYM1[name + "^-1"] = si
        assert m4mul(s, si) == m4id()
    print(f"  lambda=1: {name}: {'regular (q=1)' if reg else 'not regular with q = 1'}")
expect_not = {"chi", "chi_ac", "chi_row1", "Gm", "G5", "P", "k", "l"}
print("letters singular at 1 as expected:", expect_not.isdisjoint(SYM1), sorted(set(n for n in SYM1 if '^' not in n)))
# lambda=1 symplectic preservation and rank growth
J = [[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]]
def tr(A): return [list(r) for r in zip(*A)]
sp = all(m4mul(m4mul(tr(s), J), s) == J for s in SYM1.values())
print("lambda=1 symbols are symplectic:", sp)
names = list(SYM1); bad = 0; nI = 0
for _ in range(60):
    w = [random.choice(names) for _ in range(random.randint(2, 7))]
    s = m4id()
    for nm in w: s = m4mul(s, SYM1[nm])
    rks = [rk(w, f) for f in (16, 32, 48)]
    if s != m4id():
        if not (rks[0] < rks[1] < rks[2]): bad += 1; print("   RANK FAIL", w, rks)
    else: nI += 1
print(f"lambda=1: 60 random words, {60 - nI} with s(w) != I, rank-growth failures: {bad}")
