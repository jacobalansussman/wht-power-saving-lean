# v-meta-border, 14 Oct 2026. Counter-check of EXACTIFY.md s4 "exactify every copy ... a window of a+1 >= 2
# coefficients that EACH sub-call must then process separately ... factor c >= 2 per level, unless the unit is already
# exact at a fixed eps".  We exactify relation-using border units copy by copy with SPARSE truncated Laurent arithmetic
# (exact sub-calls), and count the paid applications actually needed (one per nonzero eps-coefficient of a paid
# move's input inside the window that can reach the output's eps^0 coefficient).
import numpy as np
ca, cb = (1 + 1j) / 2, (1 - 1j) / 2
S = np.array([[0, 1], [1, 0]], complex); D = ca * np.eye(2) + cb * S
assert np.allclose(D @ D, S)
rng = np.random.default_rng(3)
Y, N, Z = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(3)]
# Laurent-vector arithmetic: dict power -> vector (2,)
def add(*xs):
    o = {}
    for x in xs:
        for p, v in x.items(): o[p] = o.get(p, 0) + v
    return {p: v for p, v in o.items() if np.abs(v).max() > 0 or True}
def lin(Mt, x, shift=0): return {p + shift: Mt @ v for p, v in x.items()}
calls = [0]
def paid(x, prec):           # exact sub-call X = D applied coefficientwise, only to coefficients below prec
    o = {}
    for p, v in x.items():
        if p < prec: o[p] = D @ v; calls[0] += 1
    return o
v = rng.normal(size=2) + 1j * rng.normal(size=2); V = {0: v}
# (i) the toy Phi_r v = X v + eps N v + eps^-1 (X X - S) Z v   (border.exact_toy.py), M = 3 paid moves.
# Output needed mod eps^1.  X v feeds the output with coefficient eps^0 -> needs prec 1.  X(X Z v) feeds with eps^-1
# -> needs prec 2; Z v and X Z v are eps-free, so each has ONE coefficient in the window.
calls[0] = 0
u = lin(Z, V); w1 = paid(u, 2); w2 = paid(w1, 2); t = paid(V, 1)
out = add(t, lin(N, V, 1), lin(np.eye(2), add(w2, lin(-S, u)), -1))
c0 = out.get(0, 0)
print("toy Phi_r: M = 3 paid moves; coefficient-level paid calls used =", calls[0], "; |eps^0 coef - D v| =", np.abs(c0 - D @ v).max())
print("   (the eps^-1 coefficient X X Z v - S Z v is", np.abs(out.get(-1, 0)).max(), "; it never reaches eps^0, so w1, w2 can even be pruned: 1 call suffices)")
for e in [1e-1, 1e-3, 0.7]:
    Ph = D + e * N + (D @ D - S) @ Z / e
    print(f"   at fixed eps = {e}: |Phi_r(eps; D) - D| = {np.abs(Ph - D).max():.3e}  (not exact at any fixed eps != 0)")
# (ii) a two-output variant with no prunable move: out1 = X v + eps N v + eps^-1 (X X Z v - S Z v), out2 = X X Z v.
calls[0] = 0
u = lin(Z, V); w1 = paid(u, 2); w2 = paid(w1, 2); t = paid(V, 1)
out1 = add(t, lin(N, V, 1), lin(np.eye(2), add(w2, lin(-S, u)), -1)); out2 = w2
print("two-output variant: M = 3, every move needed; coefficient-level paid calls =", calls[0],
      "; errors", np.abs(out1.get(0, 0) - D @ v).max(), np.abs(out2.get(0, 0) - S @ Z @ v).max(), "-> factor 1, not >= 2")
