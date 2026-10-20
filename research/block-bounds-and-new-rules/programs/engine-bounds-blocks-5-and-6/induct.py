# pf-induction, 13 Oct 2026.  The projection recursion carried to EVERY block size (PROOF.md, Theorem 3).
# a^(m)_r = sum_s N^(m)_r(s) a^(m-1)_s  (closed-form N from nrs.py; a_0 = 0), integers, exact.
# Optional triangle improvement a_r >= a_m - (m - r) a_1 (valid for every family, PROOF.md Lemma 1(c)).
# Bases: block-2 perfect (1, 2); block-3 hook (80, 144, 224) [b56-methods]; block-4 hook k=8, f=9:
#   f_1 = 27456 (cap, exact, proved), f_2 >= 48576, f_3 >= 68320, f_4 >= 86592 (certified mod-p lower bounds,
#   b56-compute out_m4_hook8_G.txt / out_rank2_blocks.txt).
# Usage: python3 -I induct.py [M]   (default M = 60)
import sys
from fractions import Fraction as F
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nrs_closed import N_closed

def step(a, m, improve=True):
    """a = list a[0..m-1] (a[0] = 0) for block m-1  ->  list for block m."""
    b = [0] * (m + 1)
    for r in range(1, m + 1):
        tot = 0
        for s, c in N_closed(m, r).items():
            if 1 <= s <= m - 1: tot += c * a[s]
        b[r] = tot
    if improve:
        for r in range(2, m):
            b[r] = max(b[r], b[m] - (m - r) * b[1])
    return b

def run(base, M, improve=True):
    m0 = len(base); a = [0] + list(base); out = []
    for m in range(m0 + 1, M + 1):
        a = step(a, m, improve)
        out.append((m, F(a[m], a[1]), [F(a[r], a[1]) for r in range(1, m + 1)]))
    return out

if __name__ == '__main__':
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    bases = [("perfect m=2", (1, 2)), ("hook m=3 k=4", (80, 144, 224)),
             ("hook m=4 k=8", (27456, 48576, 68320, 86592))]
    for name, base in bases:
        print(f"== from {name}: profile ratios {[round(x / base[0], 4) for x in base]}")
        for improve in (True, False):
            res = run(base, M, improve)
            for m, T, prof in res:
                if m <= 8 or m % 10 == 0 or m == M:
                    tail = ", ".join(f"{float(x):.4f}" for x in prof[-4:])
                    exact = f" = {T}" if m <= 6 else ""
                    print(f"  improve={improve} m={m}: top ratio {float(T):.6f}{exact}   last four profile ratios [{tail}]")
            print(f"  improve={improve}: min over m={len(base)+1}..{M} of the top ratio = {float(min(T for _, T, _ in res)):.6f}")
