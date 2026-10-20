# an-openai: where the OpenAI unit's saving goes. Decomposition of Delta and of the dilution, per h,
# and two hypothetical ceilings: (i) side carriers removed from W (they break even, so they only dilute);
# (ii) the absolute ceiling eta <= 2*sum(dim source labels)/(W m) = 2N/(W m) with W >= 2N, i.e. eta <= 1/m.
from math import comb, log
def row(h):
    v = comb(h, 3); d = comb(h - 3, 3) + 3 * (h - 3); m = h ** 3; N = v ** 3; I = 3 * v * v; c = h + 1
    side = I * v * d; centre = I * c; bank = 2 * N; W = bank + side + centre
    L = I * c * h; Delta = 2 * N - 2 * L
    s = lambda D, WW: log(1 / (1 - D / (WW * m))) / log(m)
    return dict(h=h, m=m, bank=bank, side=side, centre=centre, W=W, gain=2 * N, loss=2 * L, Delta=Delta,
                s=s(Delta, W), s_noside=s(Delta, bank + centre), s_ceiling=s(2 * N, bank))
print("h      m     bank/W     side/W    gain 2N    loss 2L   loss/gain  saving   no-side   ceiling(1/m)")
for h in (22, 25, 30, 50, 100):
    r = row(h)
    print("%3d %7d %9.2e %9.6f %9.3e %9.3e %8.3f %9.2e %9.2e %9.2e" % (h, r["m"], r["bank"] / r["W"], r["side"] / r["W"],
          r["gain"], r["loss"], r["loss"] / r["gain"], r["s"], r["s_noside"], r["s_ceiling"]))
