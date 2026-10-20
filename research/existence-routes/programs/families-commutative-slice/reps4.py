# representations of I_4 as sums of z z^T over F_2 (4 distinct directions), and partitions of the 15 directions plus
# one extra into 4 representations (plain-cost control for U_4).
import itertools
m = 4
def cls(z):
    c = 0; b = 0
    for j in range(m):
        if z >> j & 1: c |= 1 << b
        b += 1
    for j in range(m):
        for k in range(j+1, m):
            if (z >> j & 1) and (z >> k & 1): c |= 1 << b
            b += 1
    return c
I = cls(1) ^ cls(2) ^ cls(4) ^ cls(8)
reps = [q for q in itertools.combinations(range(1, 16), 4) if cls(q[0]) ^ cls(q[1]) ^ cls(q[2]) ^ cls(q[3]) == I]
print(len(reps), "4-element representations:", reps)
for extra in range(1, 16):
    ms = list(range(1, 16)) + [extra]
    found = []
    def rec(rem, acc):
        if not rem:
            found.append(list(acc)); return
        first = rem[0]
        for r in reps:
            if r[0] != first and first not in r: continue
            if first not in r: continue
            rr = list(rem)
            ok = True
            for z in r:
                if z in rr: rr.remove(z)
                else: ok = False; break
            if ok: rec(rr, acc + [r])
    rec(sorted(ms), [])
    if found: print("extra", extra, "partitions:", found[:2], "count", len(found))
