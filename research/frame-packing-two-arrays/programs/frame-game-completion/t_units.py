"""t_units.py: lift FRAME-GAME programs (C[P] coefficients only, readouts by fpa.end_solve, exact) to engine words and
replay them on raw arrays at every f = 1..fmax (float) and, at f <= fexact, exactly modulo two primes.
Cases (block m given on the command line):
  plain1, plain2        the plain schedule, W = 1 and 2
  xchg                  r1 = y1 + s y2, r2 = y2 - s y1 (s a coordinate swap, commutes with Kbar): 2 in 2m, readout
                        K y1 = (K r1 - s K r2)/2
  o3                    the order-3 control C(2m)O3 (m >= 3): r1 = y1 + y2, r2 = y1 - s3 y2
  mix13                 r1 = y1 + y2, r2 = y1 - 2 s y2, each followed by a RANDOM free-letter path to Kbar
                        (readout with 1/3's)
  randW1 x2             random W = 1 programs q p_m T p_(m-1) ... T p_0 y (free letters with phases between moves)
Negative control: the same programs lifted with fp-build's f = 1 lift only.
usage: t_units.py m fmax fexact [fmin xmin [noneg]]"""
import sys, random, time
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
from lift import *

m = int(sys.argv[1]); fmax = int(sys.argv[2]); fexact = int(sys.argv[3])
fmin = int(sys.argv[4]) if len(sys.argv) > 4 else 1; xmin = int(sys.argv[5]) if len(sys.argv) > 5 else 1
NONEG = len(sys.argv) > 6 and sys.argv[6] == 'noneg'
E = Eng(m); S = E.S; rng = random.Random(7 * m + 1)
PR = (167772161, 469762049)

def kchain(P, r, order=None):
    for j in (order or range(m)): r = P.move(1 << j, r)
    return r
def rand_path(P, r):
    """random free letters between m moves, steered: accept only if it ends at Kbar up to a free letter"""
    for _ in range(20000):
        ps = [S.rand_P(rng) for _ in range(m + 1)]; zs = [rng.randrange(1, E.n) for _ in range(m)]
        F = ps[0]
        for t in range(m): F = S.mul(ps[t + 1], S.mul(S.T(zs[t]), F))
        if S.in_P(S.mul(S.K, S.inv(F))): break
    else: raise RuntimeError("no path")
    c = P.mix([(1, ps[0], r)])
    for t in range(m):
        c = P.move(zs[t], c); c = P.mix([(1, ps[t + 1], c)])
    return c

def cases():
    out = []
    P = FrameProg(S, 1); y = P.inp(0); kchain(P, y); out.append(("plain1", P))
    P = FrameProg(S, 2); y1, y2 = P.inp(0), P.inp(1); kchain(P, y1); kchain(P, y2); out.append(("plain2", P))
    sw = S.addr([2, 1] + [1 << k for k in range(2, m)])          # swap e0 <-> e1
    assert S.mul(sw, S.K) == S.mul(S.K, sw)
    P = FrameProg(S, 2); y1, y2 = P.inp(0), P.inp(1)
    r1 = P.mix([(1, S.I, y1), (1, sw, y2)]); r2 = P.mix([(1, S.I, y2), (-1, sw, y1)])
    kchain(P, r1); kchain(P, r2, order=list(range(m))[::-1]); out.append(("xchg", P))
    if m >= 3:
        s3 = S.addr([1 << ((k + 1) % 3) if k < 3 else 1 << k for k in range(m)])
        P = FrameProg(S, 2); y1, y2 = P.inp(0), P.inp(1)
        r1 = P.mix([(1, S.I, y1), (1, S.I, y2)]); r2 = P.mix([(1, S.I, y1), (-1, s3, y2)])
        kchain(P, r1); kchain(P, r2); out.append(("o3", P))
    P = FrameProg(S, 2); y1, y2 = P.inp(0), P.inp(1)
    r1 = P.mix([(1, S.I, y1), (1, S.I, y2)]); r2 = P.mix([(1, S.I, y1), (-2, sw, y2)])
    rand_path(P, r1); rand_path(P, r2); out.append(("mix13", P))
    for k in range(2):
        P = FrameProg(S, 1); y = P.inp(0); rand_path(P, y); out.append((f"randW1_{k}", P))
    return out

t0 = time.time()
for name, P in cases():
    st = P.solve()
    if st != 'yes':
        print(f"m={m} {name}: frame-game end solve {st} (skipped)"); continue
    nread = [len(r) for r in P.read]
    W1 = lift(E, P, 'full'); W0 = lift(E, P, 'f1')
    paid, nterms = word_stats(E, W1)
    assert paid == P.moves()
    ef = {f: float("%.1e" % run_float(E, W1, f, 100 + f)) for f in range(fmin, fmax + 1)}
    e0 = {} if NONEG else {f: float("%.1e" % run_float(E, W0, f, 100 + f)) for f in range(fmin, fmax + 1)}
    ex = {f: all(run_modp(E, W1, f, 200 + f, p) for p in PR) for f in range(xmin, fexact + 1)}
    print(f"m={m} {name}: W = {P.W}, frame moves M = {P.moves()}, readout terms {nread}; engine word: {paid} paid letters "
          f"(= M), {nterms} free terms, all free heads in B", flush=True)
    print(f"    lift (all four blocks): float errors by f {ef}")
    print(f"    exact mod {PR[0]} and {PR[1]} on random arrays, by f: {ex}")
    ex0 = {} if NONEG else {f: all(run_modp(E, W0, f, 200 + f, p) for p in PR) for f in range(xmin, min(fexact, 4) + 1)}
    print(f"    negative control, f = 1 lift only: float {e0}; exact mod p {ex0}", flush=True)
print(f"m={m}: done in {time.time() - t0:.1f} s")
