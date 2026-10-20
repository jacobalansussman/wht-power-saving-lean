# g-break g2: m = 2, f = 2. Span of M(x)M over a monoid of column-uniform zero-block steps, modulo a prime
# (i -> a square root of -1). A rank modulo a prime is a LOWER bound for the rank over Q(i). The operators
# commuting with the swap of the two columns form a space of dimension (256+16)/2 = 136 and contain K(x)K.
# So "dimension 136 reached" proves K(x)K is in the span over Q(i).
import sys, time
p = int(sys.argv[1]); mode = sys.argv[2]; cap = int(sys.argv[3])
I = next(x for x in range(2, p) if x*x % p == p-1)
d = 4
def mm(A, B): return tuple(tuple(sum(A[i][t]*B[t][j] for t in range(d)) % p for j in range(d)) for i in range(d))
def from_fn(fn):  # (M v)(x) = sum_y M[x][y] v(y)
    return tuple(tuple(fn(x, y) % p for y in range(d)) for x in range(d))
E = from_fn(lambda x, y: int(x == y))
gens = []
for z in (1, 2): gens.append(from_fn(lambda x, y, z=z: int(y == x ^ z)))                 # shifts
for G in ([2, 1], [1, 3], [3, 2]):                                                        # additive bijections: images of e1, e2
    def g(x, G=G): return (G[0] if x & 1 else 0) ^ (G[1] if x & 2 else 0)
    gens.append(from_fn(lambda x, y, g=g: int(g(x) == y)))
for z in (1, 2, 3):
    for c in (0, 1):
        gens.append(from_fn(lambda x, y, z=z, c=c: (I if (bin(x & z).count("1") + c) % 2 else 1) if x == y else 0))  # phases
fan = from_fn(lambda x, y: int(y == (x & 2)))      # fan-out: entry x reads entry (x with bit 1 cleared)
fold = from_fn(lambda x, y: int(x == (y & 2)))     # fold: entry (y with bit 1 cleared) receives the sum
if mode in ("fan", "both"): gens.append(fan)
if mode in ("fold", "both"): gens.append(fold)
def vec(M): return [M[i][j]*M[k][l] % p for i in range(d) for k in range(d) for j in range(d) for l in range(d)]
basis = {}  # pivot -> row
def add(v):
    for c in range(256):
        if v[c]:
            b = basis.get(c)
            if b is None:
                inv = pow(v[c], p-2, p); basis[c] = [x*inv % p for x in v]; return True
            f = v[c]; v = [(x - f*y) % p for x, y in zip(v, b)]
    return False
seen = {E}; frontier = [E]; add(vec(E)); t0 = time.time(); depth = 0
while frontier and len(seen) < cap and len(basis) < 136:
    new = []; depth += 1
    for A in frontier:
        for g in gens:
            B = mm(g, A)
            if B not in seen:
                seen.add(B); new.append(B); add(vec(B))
                if len(seen) >= cap or len(basis) >= 136: break
        if len(seen) >= cap or len(basis) >= 136: break
    frontier = new
# kernel K = C(x)C on the two address bits, entries with denominator 2 -> use 2K (scalar multiples do not matter)
a2, b2 = (1 + I) % p, (1 - I) % p
C2 = ((a2, b2), (b2, a2))
K = tuple(tuple(C2[(x >> 1) & 1][(y >> 1) & 1]*C2[x & 1][y & 1] % p for y in range(d)) for x in range(d))
before = len(basis); inK = not add(vec(K))
print("p=%d mode=%-5s elements=%d depth=%d span dim=%d (commutant of the swap: 136)  kernel in span mod p: %s  time %.0fs"
      % (p, mode, len(seen), depth, before, "YES" if inK else "no", time.time()-t0), flush=True)
