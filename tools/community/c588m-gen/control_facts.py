#!/usr/bin/env python3
"""m-gen: predict the Lean negative controls of C588MControls.lean with the exact mirrors (all must be FALSE):
  bad_scalar_rejected     GS.scalCheckZ badRaw (10,21,4,8,8) 0 28   (sign of the first coefficient of the first SCATTER row flipped.
                          Shea's own control, the first A gate doubled, is NOT rejected here: the first A gates of the
                          copy-free word move entrance dirt between slots and do not touch the sources 0..27)
  classic_rejected        GS.scalCheck  raw    (10,21,4,8,8) 0 28   (Sussman's classic check refuses a source write)
  duplicate_owner_rejected GX.schedK nb R badRows schedPad schedRes (row 0: entry 1 := entry 0; Shea's control)
and the positive facts they sit next to (TRUE).
usage: control_facts.py <cert.json.gz> <price-annealed.schedule.json.gz>"""
import sys, os, json, gzip, copy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zmirror as zm
from make_sched import schedK

def main(cert, sched):
    c0 = json.loads(gzip.open(cert, "rt").read()); c = zm.lean_data(c0); c["ext"] = c0["ext"]
    P = (10, 21, 4, 8, 8)
    bad = copy.deepcopy(c); k, (ng, num, den) = bad["rows"][0][0]; bad["rows"][0][0] = (k, (not ng, num, den))
    shea_style = copy.deepcopy(c); g = shea_style["A"][0]; t, (ng, num, den) = g[3][0]; g[3][0] = (t, (ng, 2 * num, den))
    s = json.loads(gzip.decompress(open(sched, "rb").read()))
    rows, pads, res = s["rows"], s["pads"], s["widths"]
    badrows = json.loads(json.dumps(rows)); badrows[0][1] = badrows[0][0]
    out = dict(first_A_gate=c0["A"][0], first_scatter_entry=c0["scat"]["table"][0][0],
               shea_style_first_A_gate_doubled_NOT_rejected=zm.evaluate(zm.scalCheck, shea_style, P, 0, 28, True),
               bad_scalar_rejected=zm.evaluate(zm.scalCheck, bad, P, 0, 28, True),
               classic_rejected=zm.evaluate(zm.scalCheck, c, P, 0, 28, False),
               positive_scal0=zm.evaluate(zm.scalCheck, c, P, 0, 28, True),
               duplicate_owner_rejected=schedK(len(pads), c["R"], badrows, pads, res),
               positive_sched=schedK(len(pads), c["R"], rows, pads, res))
    print(json.dumps(out, indent=1))
    assert not out["bad_scalar_rejected"]["value"] and "gFin" in out["bad_scalar_rejected"]["why"] and not out["classic_rejected"]["value"]
    assert out["positive_scal0"]["value"] and out["positive_sched"] and not out["duplicate_owner_rejected"]

if __name__ == "__main__":
    main(*sys.argv[1:3])
