#!/usr/bin/env python3
"""v-fp-measure: OpenAI network arithmetic, own transcription: W, M, Delta from the paper-level formulas
(W = 2N + I(vd + h + 1), M = Wm - 2N + 2L, N = v^3, I = 3v^2, L = I(h+1)h), saving ln(Wm/M)/ln m; fp-measure's closed
forms for P, D compared with its exact replay JSON (h = 6..9); ratios quoted in the NOTE."""
import json
from math import comb, log
ex = {r["h"]: r for r in json.load(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", "..", "..", "data", "packing-measurement", "out", "oa_ledger_6_7_8_9.json")))}
for h in (6, 7, 8, 9, 25, 100):
    v = comb(h, 3); d = comb(h - 3, 3) + 3 * (h - 3); m = h ** 3
    N = v ** 3; I = 3 * v * v; L = I * (h + 1) * h
    W = 2 * N + I * (v * d + h + 1); M = W * m - 2 * N + 2 * L; Dl = W * m - M
    P = v**3 * ((h-1)*(2*d+1) + 8*h) + v**3 * ((h-1)*(d+1) + h*(h-1) + 8*h + d*h*(h-1)) + v**3 * h*h*(h-1)
    D = 2*v*v*(v*d + h*(h+1) + 4*h*v) + v*v*h*(h+1)
    B = P - 2 * D - Dl
    s = log(W * m / M) / log(m) if Dl > 0 else float("nan")
    e = ex.get(h)
    cmp = ""
    if e:
        cmp = " | exact: W %s M %s P %s D %s B %s Delta %s" % (e.get("W") == W, e.get("M") == M, e.get("P") == P,
              e.get("D") == D, e.get("B") == B, e.get("Delta") == Dl)
    print("h=%d v=%d m=%d W=%.6g M=%.6g Delta=%d P/M=%.5f B/P=%.6f 2D/M=%.4g Delta/M=%.4g (P-B-2D)/Delta=%.6f maxlab v+1=%d saving=%.4g%s"
          % (h, v, m, W, M, Dl, P / M, B / P, 2 * D / M, Dl / M, (P - B - 2 * D) / Dl, v + 1, s, cmp))
