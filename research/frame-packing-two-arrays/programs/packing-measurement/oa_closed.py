#!/usr/bin/env python3
"""fp-measure: closed forms (in h) of the raw ledger of OpenAI's exchange network, read off oa_ledger.py's edge
table (edge costs as in labels_check's table; label counts as checked by oa_rows_check.py; descents from the
membership pattern measured at h = 6..9), checked against oa_ledger.py's exact counts, then evaluated at the
paper's h = 100 and at the family's best h = 25 (an-openai).  Exact integers.
  P = v^3 [ (h-1)(2d+1) + 8h ]  +  v^3 [ (h-1)(d+1) + h(h-1) + 8h + d h(h-1) ]  +  v^3 h^2 (h-1)
  D = 2 v^2 [ v d + h(h+1) + 4 h v ] + v^2 h (h+1)
  B = P - 2D - Delta,   Delta = W m - M = 2 v^3 - 2 (3 v^2)(h+1) h   (Z = I0 = 0)
usage: python3 -I oa_closed.py"""
import json, sys
from math import comb, log


def net(h):
    v = comb(h, 3)
    d = comb(h - 3, 3) + 3 * (h - 3)
    m = h ** 3
    N = v ** 3
    I = 3 * v * v
    W = 2 * N + I * (v * d + h + 1)
    L = I * (h + 1) * h
    M = W * m - 2 * N + 2 * L
    P1 = v ** 3 * ((h - 1) * (2 * d + 1) + 8 * h)
    P2 = v ** 3 * ((h - 1) * (d + 1) + h * (h - 1) + 8 * h + d * h * (h - 1))
    P3 = v ** 3 * h * h * (h - 1)
    D = 2 * v * v * (v * d + h * (h + 1) + 4 * h * v) + v * v * h * (h + 1)
    P = P1 + P2 + P3
    Delta = W * m - M
    B = P - 2 * D - Delta
    return dict(h=h, v=v, d=d, m=m, W=W, M=M, Delta=Delta, P=P, P1=P1, P2=P2, P3=P3, D=D, B=B, maxlab=v + 1)


def main():
    exact = {}
    try:
        for r in json.load(open("out/oa_ledger_6_7_8_9.json")):
            exact[r["h"]] = r
    except Exception:
        pass
    rows = []
    for h in (6, 7, 8, 9, 22, 25, 26, 100):
        r = net(h)
        e = exact.get(h)
        chk = ""
        if e:
            chk = "  exact replay: P %s, D %s, B %s" % (e["P"] == r["P"], e["D"] == r["D"], e["B"] == r["B"])
        sv = log(r["m"] * r["W"] / r["M"]) / log(r["m"]) if r["Delta"] > 0 else float("nan")
        print("h=%3d m=%d W=%.4e M=%.4e Delta=%+.4e | P/M=%.5f B/M=%.5f D/M=%.3e | B/P=%.6f | (P-B)/M=%.3e"
              " 2D/M=%.3e Delta/M=%+.3e | most labels on one move %d (c_*), next %d (targets d+2) | saving %.4g%s"
              % (h, r["m"], r["W"], r["M"], r["Delta"], r["P"] / r["M"], r["B"] / r["M"], r["D"] / r["M"],
                 r["B"] / r["P"], (r["P"] - r["B"]) / r["M"], 2 * r["D"] / r["M"], r["Delta"] / r["M"],
                 r["maxlab"], r["d"] + 2, sv, chk))
        rows.append(r)
    json.dump(rows, open("out/oa_closed.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
