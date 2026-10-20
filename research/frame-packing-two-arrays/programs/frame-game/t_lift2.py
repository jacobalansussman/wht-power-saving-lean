"""t_lift2.py: lift test with random free letters (phases included) between the moves, at block 2:
frame program c = q p3 T_z2 p2 T_z1 p1 y with q chosen so that the result is Kbar y; lifted letter by letter
(Lambda), then Pexp(K); raw arrays at f = 1..6."""
import sys, random, warnings
warnings.filterwarnings("ignore")
src = open('t_lift.py').read().split("# ---- case 1")[0].replace("m = 3; n = 1 << m", "m = 2; n = 1 << m")
exec(src)
rng = random.Random(11); found = 0
Km = np.eye(n, dtype=complex)
for j in range(m): Km = dirm(1 << j) @ Km
KP = pexp(Km)
for trial in range(400):
    ps = [S.rand_P(rng) for _ in range(3)]; zs = [rng.randrange(1, n) for _ in range(2)]
    F = ps[0]
    F = S.mul(S.T(zs[0]), F); F = S.mul(ps[1], F); F = S.mul(S.T(zs[1]), F); F = S.mul(ps[2], F)
    q = S.mul(S.K, S.inv(F))
    if not S.in_P(q): continue
    # check phases are really present
    has_phase = any(any(p[m + i] & (n - 1) for i in range(m)) for p in ps)
    errs = {}
    for f in list(range(1, 7)) + [9]:
        rn = np.random.default_rng(f)
        x = rn.normal(size=n ** f) + 1j * rn.normal(size=n ** f)
        c = Lam_P(ps[0], x, f); c = Lam_T(zs[0], c, f); c = Lam_P(ps[1], c, f); c = Lam_T(zs[1], c, f)
        c = Lam_P(ps[2], c, f); c = Lam_P(q, c, f); c = apply_free_combo(KP, c, f)
        t = tens(Km, x, f)
        errs[f] = float("%.1e" % (np.abs(c - t).max() / np.abs(t).max()))
    print(f"program {found}: phases present: {has_phase}; errors by f: {errs}", flush=True)
    found += 1
    if found == 2: break
