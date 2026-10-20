"""c3: arithmetic behind three statements of iso-lower (14 Oct 2026). Pure arithmetic, no data.

(a) Morgenstern with bound c: every gate multiplies the largest |minor| of the computed forms by at most
    max(2, c) (u+-v: 2; lam*u: |lam|). The unnormalised WHT has |det| = N^(N/2), so
        size >= (N/2) log2 N / log2 max(2, c).
    A circuit of size C * N (log2 N)^(1-s) therefore needs some |lam| >= 2^((log2 N)^s / (2C)).
    Print that requirement for our s and for the bar s = 0.1: it is vacuous at every size anyone will run,
    and only asymptotically meaningful (it beats every polylog(N) eventually).
(b) Rebalancing (Alman-Li Lemma 1.1): from depth-2 size N^(1+c), depth-2d circuits of size ~ d N^(1+c/d);
    the best d gives ~ e*c*ln N * N = (e ln2) c N log2 N. Any fixed c > 0 gives Theta(N log N).
    o(N log N) by this route needs c -> 0, i.e. depth-2 size N^(1+o(1)), which OVC/SETH forbid for H
    (Alman-Li Thm D.1).
"""
import math

print("(a) smallest admissible max|coefficient| for size C*N*(log2 N)^(1-s), C = 1")
for s in (1.6256629e-3, 1.68e-3, 0.1, 0.5):
    for logN in (20, 64, 10**3, 10**6, 10**12):
        need = (logN ** s) / 2.0            # log2 of max(2, c) must be at least this
        print(f"   s={s:<11g} log2N={logN:<14g} log2 max(2,c) >= {need:9.4f}  ->  "
              f"{'vacuous (2 already suffices)' if need <= 1 else f'c >= 2^{need:.3g}'}")

print("(b) min over d of d * N^(1 + c/d) / (N log2 N)   (c = depth-2 excess exponent)")
for c in (0.5, 0.443, 0.1, 0.01):
    for k in (20, 64, 1000):
        N = 2.0 ** k
        best = min(d * 2.0 ** (k * c / d) for d in range(1, 5 * k)) / k   # divide by log2 N (N cancels)
        print(f"   c={c:<6} log2N={k:<5} best/(N log2N) = {best:.4f}   (e ln2 c = {math.e * math.log(2) * c:.4f})")
