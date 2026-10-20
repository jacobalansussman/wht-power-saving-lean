# rl-plan: break-even head counts. Full cut (target 25, all 22,748 - t non-heads split 24+24) vs d=24 cut (target 49).
import sys; sys.path.insert(0, __file__.rsplit('/', 1)[0])
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from rescue import root, add, CALLS, W, D
HUB = D('1.6797549561255445e-3'); I = 22748
def full(t): return root(add(CALLS, 24, 2 * (I - t)), {73: W - t, 25: t})
def half(t): return root(add(add(CALLS, 24, t), 48, I - t), {73: W - t, 49: t})
lo, hi = 0, I
while hi - lo > 1:
    mid = (lo + hi) // 2
    if full(mid) > HUB: hi = mid
    else: lo = mid
print(f"full cut: break-even at t = {hi} heads ({hi/I*100:.2f} % of 22,748); at t = 11,010: {full(11010):.10e}")
for t in (1000, 5000, 11010):
    print(f"d=24 cut, t = {t:5d}: {half(t):.10e}  (+{float(half(t)-HUB):.3e})")
