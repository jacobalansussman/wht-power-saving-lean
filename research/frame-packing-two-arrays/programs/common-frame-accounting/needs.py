"""needs.py -- sl-algebra (15 Oct 2026): what a unit with saving s = 0.1 or 0.01 needs, block by block.

Saving of a unit "W in M" at block m: s = ln(mW/M)/ln m, so M/W = m^(1-s) and the move deficit per array is
delta(m, s) = m - m^(1-s).  Compared with:
  * proved lower bounds on M/W (engine game): 2 at m = 2 (Lean), 14/5 at m = 3 (Lean), 41/13 at m = 4
    (certificate-checked), 2.9447 for every m >= 4 (pf-induction, upheld);
  * the ceiling of ALL Theorem-R certificates: 169/57 at m = 3, 3.9426 at m = 4 (pf-duality D3/D4, unchecked);
  * common-frame networks whose arrays overlap their start and end frames in u dimensions on average
    (PROOF.md Theorem 1: Delta <= sum_r dim(a_r cap e_r) = u W): need u >= delta(m, s) (bridges free);
    u = 1 is every known design (lines), u = 2 is an-synth's two-sided slack.
Moves are counted one per direction (the brief's M), not at the hub's block prices.
"""
import math, json

PROVED = {2: 2.0, 3: 14 / 5, 4: 41 / 13}
TR_CEIL = {3: 169 / 57, 4: 3.9426}
OUT = {"rows": []}

def delta(m, s):
    return m - m ** (1 - s)

def first_m_line(u, s):
    """largest m with delta(m, s) <= u"""
    best = None
    for m in range(2, 5000):
        if delta(m, s) <= u:
            best = m
    return best

print(f"{'m':>4} {'M/W at 0.1':>11} {'M/W at 0.01':>12} {'proved M/W':>11} {'TR ceiling':>11} "
      f"{'0.1 status':>14} {'0.01 status':>14} {'u needed 0.1':>13} {'u needed 0.01':>14}")
for m in (2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 25, 29, 41, 73, 110):
    lb = PROVED.get(m, 2.9447 if m >= 4 else None)
    row = {"m": m, "MW_0.1": m ** 0.9, "MW_0.01": m ** 0.99, "proved_MW": lb, "TR_ceiling": TR_CEIL.get(m),
           "u_0.1": delta(m, 0.1), "u_0.01": delta(m, 0.01)}
    st1 = "excluded" if lb is not None and lb > m ** 0.9 else "open"
    st2 = "excluded" if lb is not None and lb > m ** 0.99 else "open"
    row["status_0.1"] = st1; row["status_0.01"] = st2
    OUT["rows"].append(row)
    print(f"{m:>4} {m**0.9:>11.4f} {m**0.99:>12.4f} {lb if lb is not None else float('nan'):>11.4f} "
          f"{TR_CEIL.get(m, float('nan')):>11.4f} {st1:>14} {st2:>14} {delta(m,0.1):>13.3f} {delta(m,0.01):>14.3f}")

for u in (1, 2, 4):
    for s in (0.1, 0.01):
        mm = first_m_line(u, s)
        OUT[f"largest_m_u{u}_s{s}"] = mm
        print(f"overlap u = {u} per array (bridges free): saving {s} possible only for m <= {mm}")

# W = 2 at the blocks where it would clear the bar: such units are outside the common-frame class (Corollary 2)
OUT["W2_targets"] = {"2 in 7, block 4": math.log(4 * 2 / 7) / math.log(4), "2 in 6, block 5": math.log(5 * 2 / 6) / math.log(5),
                     "2 in 5, block 3": math.log(3 * 2 / 5) / math.log(3)}
print("W = 2 targets (savings):", {k: round(v, 4) for k, v in OUT["W2_targets"].items()})
json.dump(OUT, open(__file__.replace(".py", ".json"), "w"), indent=1)
