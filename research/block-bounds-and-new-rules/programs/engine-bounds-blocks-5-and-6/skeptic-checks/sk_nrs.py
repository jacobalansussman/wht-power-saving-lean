# v-pf-induction-compute (skeptic), 13 Oct 2026. Independent code, nothing imported from induction/.
# (1) N_r(s) by brute force: for each m and r, a pair of Lagrangians at distance r built in a RANDOM symplectic
#     frame (random element of Sp applied to a standard pair), count every nonzero p by projected distance.
# (2) Compare to the closed form as stated in PROOF.md Lemma 2 (retyped here from the text).
# Vectors: ints with x = low m bits, z = high m bits; <(x|z),(x'|z')> = x.z' + x'.z.
import random, sys
def form(m, a, b):
    M = (1 << m) - 1
    return (bin(a & M & (b >> m)).count('1') + bin(b & M & (a >> m)).count('1')) & 1
def span(gens):
    S = {0}
    for g in gens:
        if g not in S: S |= {s ^ g for s in S}
    return S
def rand_sp(m, rng):
    # random symplectic map as images of the standard basis: build a random symplectic basis by
    # picking random vectors and symplectic Gram-Schmidt (only accept when full)
    n = 2 * m
    while True:
        E, F = [], []
        def red(v):
            for e, f in zip(E, F):
                if form(m, v, f): v ^= e
                if form(m, v, e): v ^= f
            return v
        ok = True
        for i in range(m):
            for _ in range(200):
                e = red(rng.randrange(1, 1 << n))
                if e: break
            else: ok = False; break
            for _ in range(200):
                f = red(rng.randrange(1, 1 << n))
                if f and form(m, e, f): break
            else: ok = False; break
            E.append(e); F.append(f)
        if ok:
            for i in range(m):
                for j in range(m):
                    assert form(m, E[i], E[j]) == 0 and form(m, F[i], F[j]) == 0 and form(m, E[i], F[j]) == (i == j)
            return E, F
def closed(m, r):
    d = {}
    for s, c in ((r, (2**(m-r) - 1) + (4**m - 2**(m+r))),
                 (r-1, 2*(2**m - 2**(m-r)) + (2**r - 1) * 2**(r-1) * 2**(m-r)),
                 (r-2, (2**r - 1) * (2**(r-1) - 1) * 2**(m-r))):
        if c: d[s] = d.get(s, 0) + c
    return d
def brute(m, r, rng):
    E, F = rand_sp(m, rng)
    L = span(E); M = span(E[:m-r] + F[m-r:])
    assert len(L & M) == 2 ** (m - r)
    cnt = {}
    for p in range(1, 1 << (2 * m)):
        A = span([v for v in L if form(m, v, p) == 0] + [p])
        B = span([v for v in M if form(m, v, p) == 0] + [p])
        assert len(A) == 2 ** m and len(B) == 2 ** m
        s = m - (len(A & B).bit_length() - 1)      # (m-1) - (dim(A cap B) - 1)
        cnt[s] = cnt.get(s, 0) + 1
    return cnt
rng = random.Random(777)
mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
bad = 0; n = 0
for m in range(1, mmax + 1):
    for r in range(1, m + 1):
        b = brute(m, r, rng); c = closed(m, r); n += 1
        if b != c: bad += 1; print("MISMATCH", m, r, b, c)
    print(f"m={m}: N_r(s) brute = closed form for r=1..{m}", flush=True)
print(f"cases {n}, mismatches {bad}")
# corrupted closed form must be rejected: move one unit of N_r(r-1) into N_r(r-2)
def closed_bad(m, r):
    d = closed(m, r)
    if r >= 2: d[r-1] -= 1; d[r-2] = d.get(r-2, 0) + 1
    return d
rej = sum(1 for m in range(2, 5) for r in range(2, m + 1) if brute(m, r, rng) != closed_bad(m, r))
print(f"corrupted closed form rejected in {rej} of {sum(1 for m in range(2,5) for r in range(2,m+1))} cases")
