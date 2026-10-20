"""T7: R8 with Z-type packing projections pi+- = (1 +- Z1)/2 and (1 +- lZ1)/2 (row-0 letters).
Kept? halves at Q? Prop S twin formula? boundedness of a seeded family at Q and at fibrewise values."""
import numpy as np, sys
sys.argv = sys.argv[:1]
exec(open('vcheck.py').read().split('SQX = np.array')[0])
rng = np.random.default_rng(777)
def rndl(d): return (rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))) / np.sqrt(2 * d)
for f in (2, 3, 4):
    L, x, r0, r1 = letters(f); n = 4 ** f; I = np.eye(n)
    Z1 = perm_mat(x, (-1.0) ** r0.sum(1))
    for nm, Pz in (("Z1", Z1), ("lZ1", L["lZ1"]), ("X1", L["X1"])):
        pp, pm = (I + Pz) / 2, (I - Pz) / 2
        res = []
        for mod in ("Q", "fib", "full"):
            if mod == "Q": D1, D2 = exactD(f, 1), exactD(f, 2)
            elif mod == "fib": D1, D2 = fib(f, rndl(2 ** f), 0), fib(f, rndl(2 ** f), 1)
            else: D1, D2 = rndl(n), rndl(n)
            u1 = D1 @ pp; u2 = D2 @ (pp + pm @ u1); u3 = D1 @ (pp @ u2)
            h = np.linalg.norm(pm @ u2) / np.sqrt(n)
            res.append(f"{mod}: |half|={h:.2f} |R8|/|half|={np.linalg.norm(pm@u2-pm@u3)/max(np.linalg.norm(pm@u2),1e-300):.1e}")
        print(f"f={f} pi=(1+-{nm})/2: " + " | ".join(res))
        # seeded family (Prop S form) and Prop S's twin formula
        dc1, dc2, dc3, dB21, dB32, dB31 = (rndl(n) for _ in range(6)); r0u = [rndl(n) for _ in range(3)]
        def family(D1, D2, s):
            u1 = D1 @ (pp + s * dc1)
            u2 = D2 @ (pp + s * dc2 + (pm + s * dB21) @ u1)
            u3 = D1 @ (s * dc3 + (pp + s * dB32) @ u2 + s * dB31 @ u1)
            return (pm @ u2 - pm @ u3) / s + r0u[0] @ u1 + r0u[1] @ u2 + r0u[2] @ u3
        D1, D2 = exactD(f, 1), exactD(f, 2)
        u1p = D1 @ (pp + pm @ dc1); u2p = D2 @ (pp + pm @ dc2 + pm @ u1p + pm @ dB21 @ pp @ u1p)
        u3p = D1 @ (pp @ u2p + pm @ dB32 @ pp @ u2p + pm @ dc3 + pm @ dB31 @ pp @ u1p)
        twin = (pm @ u2p - pm @ u3p) + r0u[0] @ pp @ u1p + r0u[1] @ pp @ u2p + r0u[2] @ pp @ u3p
        lim = [family(D1, D2, s) for s in (1e-3, 1e-5, 1e-7)]
        F1, F2 = fib(f, rndl(2 ** f), 0), fib(f, rndl(2 ** f), 1)
        print(f"      seeded at Q: |fam(1e-5)-fam(1e-7)|/|fam| = {rel(lim[1], lim[2]):.1e} (converges)  |fam - PropS twin|/|twin| = {rel(lim[2], twin):.2e}"
              f"  fibrewise |fam(1e-7)|/sqrt(n) = {np.linalg.norm(family(F1,F2,1e-7))/np.sqrt(n):.2f}")
