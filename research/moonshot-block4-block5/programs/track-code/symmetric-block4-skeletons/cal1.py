"""cal1.py -- calibration 1: H_max = B cap K B K^-1 contains a lift of every P in GL(m,2) (Levi of (L0, K L0));
tau_c fixes [K]x; tau maps the plain program to a valid program (exact span test)."""
import sys, time
from symlib import *
m = int(sys.argv[1])
S = Sym(m); G = S.G
rng = random.Random(1)
t0 = time.time()
# random GL elements
found = 0; tried = 0
Ps = []
while tried < 12:
    P = [rng.randrange(1, S.n) for _ in range(m)]
    span = {0}; ok = True
    for a in P:
        if a in span: ok = False; break
        span |= {s ^ a for s in span}
    if not ok: continue
    tried += 1
    r = S.hmax_elem(P)
    if r: found += 1; Ps.append((P, r))
print(f"m={m}: random P in GL: {found}/{tried} lift to H_max ({time.time()-t0:.1f}s)")
# tau fixes K x
for P, (c, cp) in Ps[:5]:
    kx = {(0, G.KID): Fraction(1)}
    img = S.tau(c, (0,), kx)
    print("  tau_c(Kx) == Kx:", img == kx, " c' monomial:", G.is_monomial(cp))
# tau maps the plain program to a program that still passes
W = 1
gens = plain_program(S, W)
for P, (c, cp) in Ps[:4]:
    g2 = [S.tau(c, (0,), g) for g in gens]
    # the images of the inputs must be free images of inputs
    inp_ok = all(G.is_free(gid) for (l, gid) in g2[0])
    st = final_all(S, g2, W)
    print("  tau(plain): input free:", inp_ok, "span test:", st)
# a non-H_max element breaks inputs: plain Z-phase S on qubit 0
a = [1] + [0] * (m - 1)
c = S.phase(a, [[0] * m for _ in range(m)])
cp = G.mul(S.Ki, np.conj(c).T, S.K)
print("  S-phase on qubit 0 in H_max?", G.is_monomial(cp))
