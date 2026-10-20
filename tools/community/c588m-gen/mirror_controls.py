#!/usr/bin/env python3
"""m-gen: controls of the exact mirror (zmirror.py) against facts the Lean kernel has already decided in Ryan Shea's
build (receipt of lean/complex588): Shea's certificate passes the classic checks at his parameters; his negative
control `bad_scalar_rejected` (first coefficient of the first A gate doubled) is false at batch (0, 28).
Then: the copy-free certificate FAILS the classic checks (it writes sources) and the mutable-source checks are run.
usage: mirror_controls.py <shea candidate> <c588m cert>"""
import sys, os, json, gzip, copy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zmirror as zm

def load(p):
    c0 = json.loads(gzip.open(p, "rt").read()); c = zm.lean_data(c0); c["ext"] = c0["ext"]; return c0, c

def main(shea, ours):
    out = {}
    s0, s = load(shea)
    out["shea_classic_10_21_1_8_8_batch0"] = zm.evaluate(zm.scalCheck, s, (10, 21, 1, 8, 8), 0, 28, False)
    bad = copy.deepcopy(s)
    g = bad["A"][0]; t, (ng, num, den) = g[3][0]; g[3][0] = (t, (ng, 2 * num, den))
    out["shea_bad_scalar_rejected_expect_false"] = zm.evaluate(zm.scalCheck, bad, (10, 21, 1, 8, 8), 0, 28, False)
    out["shea_classic_sw20_batch0"] = zm.evaluate(zm.scalCheck, s, (10, 20, 1, 8, 8), 0, 28, False)
    out["shea_mutable_checks_also_accept_classic_cert"] = zm.suite(s, (10, 21, 1, 8, 8), 4, (10, 8, 1, 8, 8), 80, True)["all_true"]
    o0, o = load(ours)
    out["ours_classic_expect_false"] = zm.suite(o, (10, 21, 4, 8, 8), 4, (10, 8, 4, 8, 8), 80, False)
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main(*sys.argv[1:3])
