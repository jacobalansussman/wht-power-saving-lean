#!/usr/bin/env python3
"""whatif_copies.py -- a WHAT-IF, not a result (README.md, section 7).

What would the saving be if one run of the helper circuit needed no scratch copies?

Input: a certificate as JSON (plain or .gz) with the fields h, v, R and
blocks = {"x": .., "y": .., "s": .., "c": ..}, each a map rank -> number of blocks of one invocation
(x, y: the two data arrays of a pair; s: helper arrays; c: the scratch copies).

The script builds the block list of one unit of the bridged five-stage word (m = 5h, W = 4v + R arrays): five
invocations and the idle blocks of the twins.  It finds the largest saving A/10^8 for which the whole-block rate
inequality holds at table fill 1 - 2^-13 (the convention of the eight-decimal Lean rate lemmas), once for the
list as it is and once with the class "c" deleted and everything else unchanged.  Deleting the copies raises the
deficit D = W m - (moves) by 5 * cst, cst = moves of the copies of one invocation.

Nothing here is part of a proof.  No circuit without copies is known.  Standard library only.
usage: python3 whatif_copies.py <certificate.json[.gz]> [more certificates]
"""
import sys, json, gzip
from decimal import Decimal as Dc, getcontext
from fractions import Fraction as Fr

getcontext().prec = 80


def dec(x):
    x = Fr(x)
    return Dc(x.numerator) / Dc(x.denominator)


def unit_blocks(h, v, bl, classes):
    """blocks of one unit of the bridged word: 3 forward + 2 backward invocations (same block list) + idle blocks"""
    B = {}
    for c in classes:
        for r, n in bl[c].items():
            B[int(r)] = B.get(int(r), 0) + 5 * n
    for r in (2 * (h - 1) + 4, 2 * (h - 1), h - 1, 4):       # twin 1 idle; twin 2: before, between, after
        B[r] = B.get(r, 0) + 2 * v
    return sorted(B.items())


def moment(blocks, m, W, eps, fill):
    lm = Dc(m).ln()
    z = 1 - eps
    s = sum(dec(n) * ((Dc(r).ln() - lm) * z).exp() for r, n in blocks) / dec(W)
    idle = ((Dc(m - 1).ln() - lm) * z).exp() + (-lm * z).exp()
    f = dec(fill)
    return f * s + (1 - f) * idle


def holds(blocks, m, W, S=13, scale=10 ** 8):
    """largest A with moment < 1 at saving A/scale and fill 1 - 2^-S"""
    fill = 1 - Fr(1, 2 ** S)
    lo, hi = 0, scale // 20
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if moment(blocks, m, W, Dc(mid) / scale, fill) < 1:
            lo = mid
        else:
            hi = mid
    return lo


def main(paths):
    for p in paths:
        d = json.load(gzip.open(p) if p.endswith(".gz") else open(p))
        h, v, R = d["h"], d["v"], d["R"]
        bl = d["blocks"]
        m, W = 5 * h, 4 * v + R
        cst = sum(int(r) * n for r, n in bl["c"].items())
        out = []
        for classes in ("xysc", "xys"):
            B = unit_blocks(h, v, bl, classes)
            D = W * m - sum(r * n for r, n in B)
            out.append((D, holds(B, m, W)))
        (D1, A1), (D0, A0) = out
        assert D1 == 4 * v - 5 * cst and D0 == 4 * v
        print("%s\n  h = %d, v = %d, R = %d, m = %d, W = %d, copies of one invocation %s (cst = %d moves)"
              % (p, h, v, R, m, W, dict(sorted((int(r), n) for r, n in bl["c"].items())), cst))
        print("  as it is          : D = %d, saving %d/10^8" % (D1, A1))
        print("  WHAT-IF, no copies: D = %d, saving %d/10^8  (factor %.4f)" % (D0, A0, A0 / A1))
        print("  the copies eat 5 cst / 4v = %d / %d = %.1f percent of the possible deficit"
              % (5 * cst, 4 * v, 100.0 * 5 * cst / (4 * v)))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
